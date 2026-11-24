#include "nova/physics/maxwell_cavity_energy_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace nova::physics {
namespace {
double squared_norm(const fem::EdgeVector3D& v) {
    return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}
double finite_nonnegative(double x) {
    if (!std::isfinite(x) || x<0.0)
        throw std::overflow_error("Cavity energy is not finite and nonnegative.");
    return x;
}
}

MaxwellCavityEnergy3D evaluate_maxwell_cavity_energy_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials,
    const MaxwellCavityMatrices3D& matrices,
    const MaxwellCavityMode3D& mode) {
    mesh.validate();
    const double lambda=mode.vacuum_wavenumber_squared;
    if (!(lambda>0.0) || !std::isfinite(lambda) ||
        materials.size()!=mesh.elements.size() ||
        mode.edge_field.size()!=matrices.edges.size() ||
        !std::is_sorted(matrices.edges.begin(),matrices.edges.end()) ||
        std::adjacent_find(matrices.edges.begin(),matrices.edges.end())!=matrices.edges.end())
        throw std::invalid_argument("Invalid cavity energy inputs.");
    for (const auto x:mode.edge_field)
        if (!std::isfinite(x))
            throw std::invalid_argument("Nonfinite cavity edge amplitude.");
    std::vector<bool> used(matrices.edges.size(),false);
    std::vector<bool> free(matrices.edges.size(),false);
    for (const auto e:matrices.free_edges) {
        if (e>=free.size() || free[e])
            throw std::invalid_argument("Invalid cavity free edge map.");
        free[e]=true;
    }
    for (std::size_t e=0;e<free.size();++e)
        if (!free[e] && mode.edge_field[e]!=0.0)
            throw std::invalid_argument("Nonzero constrained cavity edge.");
    MaxwellCavityEnergy3D out;
    out.cells.reserve(mesh.elements.size());
    // Degree-two symmetric tetrahedron quadrature; all weights positive.
    constexpr double a=0.5854101966249685, b=0.1381966011250105;
    constexpr std::array<std::array<double,3>,4> quadrature{{
        {b,b,b},{a,b,b},{b,a,b},{b,b,a}}};
    constexpr std::array<std::array<double,3>,4> vertices{{
        {0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        const auto& material=materials[c];
        const double eps=material.relative_permittivity.real();
        const double mu=material.relative_permeability;
        if (!(eps>0.0) || !std::isfinite(eps) ||
            material.relative_permittivity.imag()!=0.0 ||
            !(mu>0.0) || !std::isfinite(mu))
            throw std::invalid_argument("Cavity energy requires lossless positive scalar materials.");
        const auto cell=fem::make_nedelec_tetrahedron_3d(mesh,c);
        std::array<double,6> coefficients{};
        fem::EdgeVector3D curl{};
        for (std::size_t i=0;i<6;++i) {
            const auto it=std::lower_bound(matrices.edges.begin(),matrices.edges.end(),cell.edges[i]);
            if (it==matrices.edges.end() || *it!=cell.edges[i])
                throw std::invalid_argument("Missing cavity mesh edge.");
            const auto e=static_cast<std::size_t>(it-matrices.edges.begin());
            used[e]=true;
            coefficients[i]=mode.edge_field[e];
            for (std::size_t d=0;d<3;++d) curl[d]+=coefficients[i]*cell.curls[i][d];
        }
        const auto electric_weight=[&](const std::array<double,3>& p) {
            const auto basis=fem::evaluate_nedelec_basis_3d(cell,p[0],p[1],p[2]);
            fem::EdgeVector3D field{};
            for (std::size_t i=0;i<6;++i)
                for (std::size_t d=0;d<3;++d)
                    field[d]+=coefficients[i]*basis[i][d];
            return finite_nonnegative(eps*squared_norm(field));
        };
        MaxwellCavityCellEnergy3D energy;
        for (const auto& p:quadrature)
            energy.electric_integral+=electric_weight(p)*(cell.volume/4.0);
        energy.electric_integral=finite_nonnegative(energy.electric_integral);
        energy.magnetic_integral=finite_nonnegative(
            ((squared_norm(curl)/mu)/lambda)*cell.volume);
        for (const auto& p:vertices)
            out.peak_electric_weight=std::max(out.peak_electric_weight,electric_weight(p));
        out.electric_integral+=energy.electric_integral;
        out.magnetic_integral+=energy.magnetic_integral;
        out.cells.push_back(energy);
    }
    if (std::find(used.begin(),used.end(),false)!=used.end())
        throw std::invalid_argument("Unexpected cavity mesh edge.");
    const double total=finite_nonnegative(out.electric_integral+out.magnetic_integral);
    if (!(out.electric_integral>0.0) || !(out.peak_electric_weight>0.0))
        throw std::invalid_argument("Zero or underflowed cavity electric energy.");
    out.relative_energy_imbalance=std::abs(out.electric_integral-out.magnetic_integral)/total;
    out.electric_mode_volume_m3=finite_nonnegative(out.electric_integral/out.peak_electric_weight);
    for (auto& cell:out.cells) {
        cell.electric_participation=cell.electric_integral/out.electric_integral;
        cell.total_participation=(cell.electric_integral+cell.magnetic_integral)/total;
    }
    return out;
}
} // namespace nova::physics
