#include "nova/physics/maxwell_cavity_3d.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
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
} // namespace nova::physics
