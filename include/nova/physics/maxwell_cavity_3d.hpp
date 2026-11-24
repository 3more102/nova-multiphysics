#pragma once
#include "nova/physics/maxwell_3d.hpp"
#include "nova/linalg/symmetric_eigen.hpp"

namespace nova::physics {

struct MaxwellCavityMatrices3D {
    std::vector<fem::Edge3D> edges;
    std::vector<std::size_t> free_edges;
    linalg::CsrMatrix stiffness,mass;
};

// PEC on the entire boundary, lossless positive scalar epsilon_r and mu_r.
[[nodiscard]] MaxwellCavityMatrices3D assemble_maxwell_cavity_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellMaterial3D> materials);

struct MaxwellCavityMode3D {
    double vacuum_wavenumber_squared{};
    double frequency_hz{};
    std::vector<double> edge_field; // full edge ordering, zero on PEC
    double relative_residual{};
};
struct MaxwellCavityResult3D {
    std::vector<MaxwellCavityMode3D> modes;
    std::size_t nullspace_dimension{};
    double nullspace_cutoff{};
    std::size_t sweeps{};
    bool converged{};
};

// Numerically excludes zero curl-gradient eigenvalues. This is a small
// reference solver, not a scalable eigensolver or general spurious-mode proof.
[[nodiscard]] MaxwellCavityResult3D solve_maxwell_cavity_3d(
    const MaxwellCavityMatrices3D& matrices,std::size_t max_modes,
    double relative_nullspace_tolerance=1e-10,
    const linalg::SymmetricEigenOptions& options={});

} // namespace nova::physics
