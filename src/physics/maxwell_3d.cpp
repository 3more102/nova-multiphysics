#include "nova/physics/maxwell_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

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
    std::set<fem::Edge3D> pec;
    if (pec_on_exterior)
        for (const auto& face : fem::extract_boundary_faces_3d(mesh))
            for (std::size_t i=0;i<3;++i)
                for (std::size_t j=i+1;j<3;++j)
                    pec.insert({std::min(face[i],face[j]),std::max(face[i],face[j])});
    MaxwellSystem3D out;
    const auto invalid=std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> reduced(map.size(),invalid);
    for (const auto& [edge,id] : map) {
        out.edges.push_back(edge);
        if (!pec.contains(edge)) {
            reduced[id]=out.free_edges.size();
            out.free_edges.push_back(id);
        }
    }
    const auto n=out.free_edges.size();
    if (n==0 || n>invalid/2)
        throw std::invalid_argument("Maxwell system has no free edges or is too large.");
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

} // namespace nova::physics
