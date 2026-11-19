#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

// Experimental P1 tetrahedral, time-independent-coefficient transient heat model:
//   C(x) dT/dt - div(K(x) grad(T)) = Q(x).
// Reuses the stationary 3D material, volumetric source and boundary APIs.
// All boundary functions and loads are fixed in time for this first 3D release.
struct TransientHeat3DProblem {
    ScalarDiffusion3DProblem diffusion;
    // Exactly one representation: positive spatial C(x) or one positive C per
    // tetrahedron, in mesh.elements ordering. C has units J/(m^3 K).
    ScalarField3D volumetric_capacity;
    std::vector<double> cell_volumetric_capacity;
    // Exact initial nodal field; Dirichlet values are imposed at t=start_time.
    std::vector<double> initial_temperature;
};

enum class TransientHeat3DScheme {
    BackwardEuler,
    CrankNicolson
};

struct TransientHeatStepDiagnostics3D;

struct TransientHeat3DOptions {
    TransientHeat3DScheme scheme{TransientHeat3DScheme::BackwardEuler};
    double time_step{0.1};
    double start_time{0.0};
    std::size_t steps{1U};
    linalg::SolverOptions solver_options{};
    // True preserves legacy per-step full temperature vectors in linear_solvers.
    // False keeps solver metadata and final temperature, but drops these vectors.
    bool retain_step_solution_vectors{true};
    // Called after each successful step with the physical state and energy budget.
    // The temperature reference is valid only during the callback; copy if needed.
    // Exceptions from the callback propagate and abort the solve.
    std::function<void(std::size_t, const std::vector<double>&,
                       const TransientHeatStepDiagnostics3D&)> on_step{};
};

// All energies are integrated over a single time step in joules.
// Neumann energy is negative for imposed outward heat; Robin energy includes
// transfer to/from ambient. Dirichlet energy is the discrete reaction at
// prescribed nodes, NOT an equilibrated boundary-face flux.
struct TransientHeatStepDiagnostics3D {
    double time{};
    double heat_content{};    // Integral C*T dV (J)
    double source_energy{};   // Volumetric source input (J)
    double neumann_energy{};  // Prescribed surface flux input (J)
    double robin_energy{};    // Net convection input (J)
    double dirichlet_energy{};// Reaction at temperature-prescribed DOFs (J)
    double balance_residual{};// Discrete heat-budget closure (J)
};

struct TransientHeat3DResult {
    std::vector<double> temperature;
    double final_time{};
    std::vector<linalg::SolverResult> linear_solvers;
    std::vector<TransientHeatStepDiagnostics3D> diagnostics;
};

// Consistent P1 thermal mass and existing validated 3D conductivity integration.
// One time-independent system matrix is reused across all steps; CG receives a
// warm-start state. Purely insulated meshes are permitted because positive
// thermal mass removes the steady-state constant-mode singularity.
// No time-varying coefficients, phase change or adaptive 3D time steps yet.
[[nodiscard]] TransientHeat3DResult solve_transient_heat_3d(
    const TransientHeat3DProblem& problem,
    const TransientHeat3DOptions& options = {}
);

} // namespace nova::fem
