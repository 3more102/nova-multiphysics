#include "nova/physics/maxwell_cavity_3d.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace nova::physics {
MaxwellCavityResult3D solve_maxwell_cavity_3d(
    const MaxwellCavityMatrices3D& matrices,std::size_t max_modes,
    double relative_nullspace_tolerance,const linalg::SymmetricEigenOptions& options) {
    if (max_modes==0 || !(relative_nullspace_tolerance>0.0)
        || relative_nullspace_tolerance>=1.0 || !std::isfinite(relative_nullspace_tolerance)
        || matrices.free_edges.size()!=matrices.stiffness.rows())
        throw std::invalid_argument("Invalid cavity options or dimensions.");
    for (std::size_t i=0;i<matrices.free_edges.size();++i)
        if (matrices.free_edges[i]>=matrices.edges.size() ||
            (i>0 && matrices.free_edges[i]<=matrices.free_edges[i-1]))
            throw std::invalid_argument("Invalid cavity edge ordering.");
    const auto eig=linalg::solve_symmetric_generalized_eigen(
        matrices.stiffness,matrices.mass,options);
    MaxwellCavityResult3D out;
    out.converged=eig.converged; out.sweeps=eig.sweeps;
    double scale=0.0;
    for (const auto& p:eig.eigenpairs) scale=std::max(scale,std::abs(p.eigenvalue));
    out.nullspace_cutoff=scale*relative_nullspace_tolerance;
    constexpr double c=299792458.0;
    for (const auto& p:eig.eigenpairs) {
        if (p.eigenvalue < -out.nullspace_cutoff)
            throw std::runtime_error("Cavity curl operator has negative eigenvalue.");
        if (std::abs(p.eigenvalue)<=out.nullspace_cutoff) {
            ++out.nullspace_dimension; continue;
        }
        if (out.modes.size()==max_modes) continue;
        MaxwellCavityMode3D mode;
        mode.vacuum_wavenumber_squared=p.eigenvalue;
        mode.frequency_hz=c/(2*std::numbers::pi)*std::sqrt(p.eigenvalue);
        if (!std::isfinite(mode.frequency_hz)) throw std::overflow_error("Cavity frequency overflow.");
        mode.relative_residual=p.relative_residual;
        mode.edge_field.resize(matrices.edges.size());
        for (std::size_t i=0;i<matrices.free_edges.size();++i)
            mode.edge_field[matrices.free_edges[i]]=p.eigenvector[i];
        out.modes.push_back(std::move(mode));
    }
    return out;
}

MaxwellCavityParticipation3D maxwell_cavity_participation_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    const MaxwellCavityMatrices3D& matrices,
    std::span<const double> edge_field) {
    mesh.validate();
    std::map<fem::Edge3D,std::size_t> map;
    for (const auto& element:mesh.elements)
        for (const auto& [i,j]:fem::tetrahedron_local_edges)
            map.try_emplace({std::min(element[i],element[j]),std::max(element[i],element[j])},0);
    std::set<fem::Edge3D> pec;
    for (const auto& face:fem::extract_boundary_faces_3d(mesh))
        for (std::size_t i=0;i<3;++i)
            for (std::size_t j=i+1;j<3;++j)
                pec.insert({std::min(face[i],face[j]),std::max(face[i],face[j])});
    std::vector<fem::Edge3D> edges;
    std::vector<std::size_t> free_edges;
    for (auto& [edge,index]:map) {
        index=edges.size();
        edges.push_back(edge);
        if (!pec.contains(edge)) free_edges.push_back(index);
    }
    if (materials.size()!=mesh.elements.size() ||
        edges!=matrices.edges || free_edges!=matrices.free_edges ||
        edge_field.size()!=edges.size())
        throw std::invalid_argument("Cavity participation topology or size mismatch.");
    for (std::size_t i=0;i<edge_field.size();++i) {
        if (!std::isfinite(edge_field[i]) ||
            (!std::binary_search(free_edges.begin(),free_edges.end(),i)
             && edge_field[i]!=0.0))
            throw std::invalid_argument("Cavity participation requires finite PEC coefficients.");
    }
    MaxwellCavityParticipation3D out;
    constexpr double a=0.5854101966249685,b=0.1381966011250105;
    constexpr std::array<std::array<double,3>,4> points{{
        {b,b,b},{a,b,b},{b,a,b},{b,b,a}
    }};
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        const auto& m=materials[c];
        if (!std::isfinite(m.relative_permittivity.real()) || m.relative_permittivity.imag()!=0.0 ||
            !(m.relative_permittivity.real()>0.0) ||
            !(m.relative_permeability>0.0) || !std::isfinite(m.relative_permeability))
            throw std::invalid_argument("Cavity participation requires lossless positive materials.");
        const auto cell=fem::make_nedelec_tetrahedron_3d(mesh,c);
        std::array<double,6> coeff{};
        fem::EdgeVector3D curl{};
        for (std::size_t e=0;e<6;++e) {
            coeff[e]=edge_field[map.at(cell.edges[e])];
            for (std::size_t d=0;d<3;++d) curl[d]+=coeff[e]*cell.curls[e][d];
        }
        double mass=0.0,stiffness=0.0;
        for (const auto& p:points) {
            const auto basis=fem::evaluate_nedelec_basis_3d(cell,p[0],p[1],p[2]);
            fem::EdgeVector3D field{};
            for (std::size_t e=0;e<6;++e)
                for (std::size_t d=0;d<3;++d) field[d]+=coeff[e]*basis[e][d];
            for (double v:field) mass+=v*v;
        }
        mass*=cell.volume/4.0*m.relative_permittivity.real();
        for (double v:curl) stiffness+=v*v;
        stiffness*=cell.volume/m.relative_permeability;
        if (!std::isfinite(mass) || !std::isfinite(stiffness))
            throw std::overflow_error("Cavity participation energy overflow.");
        out.electric_fraction.push_back(mass);
        out.magnetic_fraction.push_back(stiffness);
        out.mass_norm_squared+=mass;
        out.curl_energy+=stiffness;
    }
    if (!std::isfinite(out.mass_norm_squared) || !std::isfinite(out.curl_energy))
        throw std::overflow_error("Cavity participation total overflow.");
    if (!(out.mass_norm_squared>0.0) || !(out.curl_energy>0.0))
        throw std::invalid_argument("Cavity participation requires nonzero electric and curl energy.");
    for (auto& value:out.electric_fraction) value/=out.mass_norm_squared;
    for (auto& value:out.magnetic_fraction) value/=out.curl_energy;
    return out;
}

} // namespace nova::physics
