#include "nova/physics/maxwell_3d.hpp"
#include "nova/physics/maxwell_cavity_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace nova::physics {
namespace {
using Complex=std::complex<double>;
bool finite(Complex z) {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}
using EdgeMap=std::map<fem::Edge3D,std::size_t>;
EdgeMap edge_map(const fem::TetrahedralMesh3D& mesh) {
    EdgeMap map;
    for (const auto& cell : mesh.elements)
        for (const auto& [i,j] : fem::tetrahedron_local_edges)
            map.try_emplace({std::min(cell[i],cell[j]),std::max(cell[i],cell[j])},0);
    std::size_t index=0;
    for (auto& [edge,id] : map) {
        (void)edge;
        id=index++;
    }
    return map;
}
MaxwellSystem3D make_topology(const fem::TetrahedralMesh3D& mesh,
    const EdgeMap& map,bool pec_on_exterior) {
    std::set<fem::Edge3D> pec;
    if (pec_on_exterior)
        for (const auto& face : fem::extract_boundary_faces_3d(mesh))
            for (std::size_t i=0;i<3;++i)
                for (std::size_t j=i+1;j<3;++j)
                    pec.insert({std::min(face[i],face[j]),std::max(face[i],face[j])});
    MaxwellSystem3D out;
    const auto invalid=std::numeric_limits<std::size_t>::max();
    for (const auto& [edge,id] : map) {
        out.edges.push_back(edge);
        if (!pec.contains(edge)) {
            out.free_edges.push_back(id);
        }
    }
    const auto n=out.free_edges.size();
    if (n==0 || n>invalid/2)
        throw std::invalid_argument("Maxwell system has no free edges or is too large.");
    return out;
}
std::vector<std::size_t> make_reduced(const MaxwellSystem3D& topology) {
    std::vector<std::size_t> reduced(topology.edges.size(),
        std::numeric_limits<std::size_t>::max());
    for (std::size_t i=0;i<topology.free_edges.size();++i)
        reduced[topology.free_edges[i]]=i;
    return reduced;
}
}

MaxwellSystem3D assemble_maxwell_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    double vacuum_wavenumber, bool pec_on_exterior) {
    mesh.validate();
    const double k2=vacuum_wavenumber*vacuum_wavenumber;
    if (!(vacuum_wavenumber>0.0) || !(k2>0.0) || !std::isfinite(k2)
        || materials.size()!=mesh.elements.size())
        throw std::invalid_argument("Invalid Maxwell frequency or material count.");
    for (const auto& m : materials)
        if (!finite(m.relative_permittivity) || m.relative_permittivity.imag()<0.0
            || !(m.relative_permeability>0.0) || !std::isfinite(m.relative_permeability))
            throw std::invalid_argument("Invalid passive Maxwell cell material.");
    const auto map=edge_map(mesh);
    auto out=make_topology(mesh,map,pec_on_exterior);
    const auto n=out.free_edges.size();
    const auto invalid=std::numeric_limits<std::size_t>::max();
    const auto reduced=make_reduced(out);
    std::vector<linalg::Triplet> entries;
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        const auto cell=fem::make_nedelec_tetrahedron_3d(mesh,c);
        const auto local=fem::nedelec_element_matrices_3d(cell);
        const auto& m=materials[c];
        for (std::size_t e=0;e<6;++e) {
            const auto i=reduced[map.at(cell.edges[e])];
            if (i==invalid) continue;
            for (std::size_t f=0;f<6;++f) {
                const auto j=reduced[map.at(cell.edges[f])];
                if (j==invalid) continue;
                const Complex a=local.curl_curl[e][f]/m.relative_permeability
                              -k2*local.mass[e][f]*m.relative_permittivity;
                if (!finite(a)) throw std::overflow_error("Maxwell matrix overflow.");
                if (a.real()!=0.0) {
                    entries.push_back({i,j,a.real()});
                    entries.push_back({i+n,j+n,a.real()});
                }
                if (a.imag()!=0.0) {
                    entries.push_back({i,j+n,-a.imag()});
                    entries.push_back({i+n,j,a.imag()});
                }
            }
        }
    }
    out.real_block=linalg::CsrMatrix::from_triplets(2*n,2*n,entries);
    return out;
}

std::vector<Complex> assemble_maxwell_cell_load_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const ComplexVector3D> cell_sources,
    const MaxwellSystem3D& system) {
    mesh.validate();
    const auto map=edge_map(mesh);
    if (cell_sources.size()!=mesh.elements.size() || system.edges.size()!=map.size())
        throw std::invalid_argument("Maxwell source or edge count mismatch.");
    for (const auto& [edge,id] : map)
        if (system.edges[id]!=edge) throw std::invalid_argument("Maxwell edge ordering mismatch.");
    std::vector<Complex> load(map.size());
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        for (const auto z : cell_sources[c])
            if (!finite(z)) throw std::invalid_argument("Nonfinite Maxwell source.");
        const auto cell=fem::make_nedelec_tetrahedron_3d(mesh,c);
        // Centroid rule is exact for the affine basis times constant source.
        const auto basis=fem::evaluate_nedelec_basis_3d(cell,0.25,0.25,0.25);
        for (std::size_t e=0;e<6;++e) {
            auto& value=load[map.at(cell.edges[e])];
            for (std::size_t d=0;d<3;++d)
                value+=cell.volume*basis[e][d]*cell_sources[c][d];
            if (!finite(value)) throw std::overflow_error("Maxwell load overflow.");
        }
    }
    return load;
}

MaxwellSolution3D solve_maxwell_3d(
    const MaxwellSystem3D& system, std::span<const Complex> edge_load,
    const linalg::GmresOptions& options) {
    const auto n=system.free_edges.size();
    if (n==0 || n>std::numeric_limits<std::size_t>::max()/2 ||
        system.real_block.rows()!=2*n || system.real_block.columns()!=2*n ||
        edge_load.size()!=system.edges.size())
        throw std::invalid_argument("Invalid Maxwell system or load dimensions.");
    for (auto z : edge_load)
        if (!finite(z)) throw std::invalid_argument("Nonfinite Maxwell load.");
    std::vector<double> rhs(2*n);
    for (std::size_t i=0;i<n;++i) {
        const auto index=system.free_edges[i];
        if (index>=edge_load.size() || (i>0 && index<=system.free_edges[i-1]))
            throw std::invalid_argument("Invalid Maxwell free edge ordering.");
        rhs[i]=edge_load[index].real();
        rhs[i+n]=edge_load[index].imag();
    }
    MaxwellSolution3D out;
    out.diagnostics=linalg::solve_gmres(system.real_block,rhs,options);
    if (out.diagnostics.solution.size()!=2*n)
        throw std::runtime_error("Unexpected GMRES solution dimension.");
    out.edge_field.resize(system.edges.size());
    for (std::size_t i=0;i<n;++i)
        out.edge_field[system.free_edges[i]]={out.diagnostics.solution[i],out.diagnostics.solution[i+n]};
    return out;
}


namespace {
// SI vacuum constants. c is exact; epsilon0 uses the 2018 CODATA value.
// mu0 follows from epsilon0*mu0*c^2=1.
constexpr double light_speed=299792458.0;
constexpr double vacuum_epsilon=8.8541878128e-12;
constexpr double vacuum_mu=1.0/(vacuum_epsilon*light_speed*light_speed);

double checked_omega(const MaxwellCellMaterial3D& m,double k) {
    const double omega=k*light_speed;
    if (!(omega>0.0) || !std::isfinite(omega) || !finite(m.relative_permittivity)
        || m.relative_permittivity.imag()<0.0 || !(m.relative_permeability>0.0)
        || !std::isfinite(m.relative_permeability))
        throw std::invalid_argument("Invalid Maxwell field material or frequency.");
    return omega;
}
}

MaxwellFieldSample3D evaluate_maxwell_field_3d(
    const fem::NedelecTetrahedron3D& cell,
    const std::array<Complex,6>& coefficients,
    const MaxwellCellMaterial3D& material, double vacuum_wavenumber,
    double r,double s,double t) {
    const double omega=checked_omega(material,vacuum_wavenumber);
    const double magnetic_scale=omega*vacuum_mu*material.relative_permeability;
    if (!(magnetic_scale>0.0) || !std::isfinite(magnetic_scale))
        throw std::overflow_error("Invalid magnetic field scale.");
    const auto basis=fem::evaluate_nedelec_basis_3d(cell,r,s,t);
    MaxwellFieldSample3D out;
    ComplexVector3D curl{};
    for (std::size_t e=0;e<6;++e) {
        if (!finite(coefficients[e])) throw std::invalid_argument("Nonfinite edge coefficient.");
        for (std::size_t d=0;d<3;++d) {
            out.electric_field[d]+=coefficients[e]*basis[e][d];
            curl[d]+=coefficients[e]*cell.curls[e][d];
        }
    }
    for (std::size_t d=0;d<3;++d) {
        out.magnetic_field[d]=curl[d]/Complex(0.0,magnetic_scale);
        if (!finite(out.electric_field[d]) || !finite(out.magnetic_field[d]))
            throw std::overflow_error("Maxwell field reconstruction overflow.");
    }
    for (std::size_t d=0;d<3;++d) {
        const auto j=(d+1)%3, k=(d+2)%3;
        out.time_averaged_poynting[d]=0.5*std::real(
            out.electric_field[j]*std::conj(out.magnetic_field[k])
           -out.electric_field[k]*std::conj(out.magnetic_field[j]));
        if (!std::isfinite(out.time_averaged_poynting[d]))
            throw std::overflow_error("Maxwell Poynting flux overflow.");
    }
    if (material.relative_permittivity.imag()>0.0) {
        double intensity=0.0;
        for (auto z : out.electric_field) intensity+=std::norm(z);
        out.absorbed_power_density=0.5*omega*vacuum_epsilon
            *material.relative_permittivity.imag()*intensity;
        if (!std::isfinite(out.absorbed_power_density))
            throw std::overflow_error("Maxwell dielectric loss overflow.");
    }
    return out;
}

std::vector<double> integrate_maxwell_absorption_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    double vacuum_wavenumber,
    const MaxwellSystem3D& system, std::span<const Complex> edge_field) {
    mesh.validate();
    const auto map=edge_map(mesh);
    if (materials.size()!=mesh.elements.size() || map.size()!=system.edges.size()
        || edge_field.size()!=map.size())
        throw std::invalid_argument("Maxwell absorption dimensions mismatch.");
    for (const auto& [edge,id] : map)
        if (system.edges[id]!=edge)
            throw std::invalid_argument("Maxwell absorption edge ordering mismatch.");
    // Symmetric four-point tetrahedral rule, degree two, positive weights.
    constexpr double a=0.5854101966249685, b=0.1381966011250105;
    constexpr std::array<std::array<double,3>,4> points{{
        {b,b,b},{a,b,b},{b,a,b},{b,b,a}
    }};
    std::vector<double> result(mesh.elements.size());
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        const auto cell=fem::make_nedelec_tetrahedron_3d(mesh,c);
        std::array<Complex,6> coefficients{};
        for (std::size_t e=0;e<6;++e)
            coefficients[e]=edge_field[map.at(cell.edges[e])];
        for (const auto& p : points) {
            const auto sample=evaluate_maxwell_field_3d(
                cell,coefficients,materials[c],vacuum_wavenumber,p[0],p[1],p[2]);
            result[c]+=cell.volume/4.0*sample.absorbed_power_density;
        }
        if (!std::isfinite(result[c]))
            throw std::overflow_error("Maxwell integrated loss overflow.");
    }
    return result;
}


MaxwellCavityMatrices3D assemble_maxwell_cavity_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials) {
    mesh.validate();
    if (materials.size()!=mesh.elements.size())
        throw std::invalid_argument("Cavity material count mismatch.");
    for (const auto& m : materials)
        if (!finite(m.relative_permittivity) || m.relative_permittivity.imag()!=0.0
            || !(m.relative_permittivity.real()>0.0)
            || !(m.relative_permeability>0.0) || !std::isfinite(m.relative_permeability))
            throw std::invalid_argument("Cavity requires lossless positive isotropic materials.");
    const auto map=edge_map(mesh);
    auto topology=make_topology(mesh,map,true);
    const auto reduced=make_reduced(topology);
    const auto invalid=std::numeric_limits<std::size_t>::max();
    const auto n=topology.free_edges.size();
    std::vector<linalg::Triplet> stiffness,mass;
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        const auto cell=fem::make_nedelec_tetrahedron_3d(mesh,c);
        const auto local=fem::nedelec_element_matrices_3d(cell);
        for (std::size_t e=0;e<6;++e) {
            const auto i=reduced[map.at(cell.edges[e])];
            if (i==invalid) continue;
            for (std::size_t f=0;f<6;++f) {
                const auto j=reduced[map.at(cell.edges[f])];
                if (j==invalid) continue;
                const double k=local.curl_curl[e][f]/materials[c].relative_permeability;
                const double m=local.mass[e][f]*materials[c].relative_permittivity.real();
                if (!std::isfinite(k) || !std::isfinite(m))
                    throw std::overflow_error("Cavity matrix overflow.");
                stiffness.push_back({i,j,k});
                mass.push_back({i,j,m});
            }
        }
    }
    return {std::move(topology.edges),std::move(topology.free_edges),
        linalg::CsrMatrix::from_triplets(n,n,stiffness),
        linalg::CsrMatrix::from_triplets(n,n,mass)};
}

} // namespace nova::physics
