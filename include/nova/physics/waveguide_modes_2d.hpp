#pragma once

#include "nova/fem/mesh2d.hpp"
#include "nova/linalg/symmetric_eigen.hpp"

#include <cstddef>
#include <vector>

namespace nova::physics {

// Scalar z-invariant reductions. The names identify the solved field:
// TE_Ez solves for E_z; TM_Hz solves for H_z.
enum class WaveguidePolarization2D {
    TE_Ez,
    TM_Hz
};

struct WaveguideMaterial2D {
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
};

struct WaveguideModeProblem2D {
    fem::TriangleMesh2D mesh;
    std::vector<WaveguideMaterial2D> materials; // one lossless material per triangle
    double vacuum_wavelength_m{};
    double cladding_index{};
    WaveguidePolarization2D polarization{WaveguidePolarization2D::TE_Ez};

    // Homogeneous scalar-field Dirichlet constraints. An empty list constrains
    // every exterior mesh vertex. Unconstrained exterior edges have natural
    // homogeneous Neumann conditions; this permits symmetry/truncation planes.
    std::vector<std::size_t> dirichlet_boundary_nodes;
};

struct WaveguideMode2D {
    double effective_index{};
    double propagation_constant_per_m{};
    double propagation_constant_squared_per_m2{};
    double core_confinement_factor{}; // weighted scalar energy in cells above the cladding index
    std::vector<double> scalar_field; // full mesh-vertex ordering; zero on constrained nodes
    double relative_residual{};
};

struct WaveguideModesResult2D {
    std::vector<WaveguideMode2D> modes; // descending effective-index order
    std::size_t constrained_nodes{};
    std::size_t free_nodes{};
    std::size_t jacobi_sweeps{};
    bool converged{};
};

// P1 triangular FEM for lossless, isotropic, piecewise-constant materials.
// The bounded dense generalized eigensolver currently limits the free space
// to 512 DOFs. Light-line candidates with cladding_index < n_eff < max(n_material)
// are returned with a core-confinement diagnostic; finite-window continuum modes
// may also appear. This is a scalar reference solver, not full-vector Maxwell FEM.
[[nodiscard]] WaveguideModesResult2D solve_waveguide_modes_2d(
    const WaveguideModeProblem2D& problem,
    std::size_t max_modes = 6,
    const linalg::SymmetricEigenOptions& eigen_options = {});

} // namespace nova::physics
