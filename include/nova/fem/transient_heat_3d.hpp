#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <functional>
#include <span>
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
    // Defaults preserve the full per-step solver/diagnostic history.
    // Disable for large runs to bound solver working memory by mesh/DOFs, not time steps.
    bool store_step_history{true};
    // Called for every completed step, including when history storage is off.
    // The temperature span and diagnostics reference are only valid during
    // the callback. Exceptions propagate and abort the solve.
    std::function<void(std::size_t step, double time,
                       std::span<const double> temperature,
                       const TransientHeatStepDiagnostics3D& diagnostics)> on_step{};
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
    // Empty when options.store_step_history is false.
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
