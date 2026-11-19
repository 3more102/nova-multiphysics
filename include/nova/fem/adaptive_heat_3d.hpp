#pragma once

#include "nova/fem/transient_heat_3d.hpp"

#include <cstddef>
#include <vector>

namespace nova::fem {

// Fixed-mesh adaptive time stepping for linear P1 tetrahedral heat transfer.
// Each trial compares one full step with two half steps from the same state.
// The accepted state is the two-half-step solution.
struct AdaptiveHeat3DOptions {
    TransientTimeScheme scheme{TransientTimeScheme::CrankNicolson};
    double start_time{0.0};
    double end_time{1.0};
    double initial_step{0.1};
    double min_step{1.0e-6};
    double max_step{0.25};
    double absolute_tolerance{1.0e-6};
    double relative_tolerance{1.0e-4};
    std::size_t max_attempts{10000U};
    linalg::SolverOptions solver_options{};
};

struct AdaptiveHeatAttempt3D {
    double start_time{};
    double attempted_step{};
    // RMS over all P1 vertex DOFs of the step-doubling estimate:
    // |T_two_half - T_full|/(2^p - 1), scaled by atol + rtol*max(|T_old|,|T_two_half|).
    // A value <= 1 permits acceptance.
    double normalized_error{};
    bool accepted{};
};

struct AdaptiveHeat3DResult {
    std::vector<double> temperature;
    double final_time{};
    std::size_t accepted_steps{};
    std::size_t rejected_steps{};
    // IC(0) cache setup counts, including changes from adaptive dt or Robin
    // coefficients. Both are zero when IC(0) is not selected.
    std::size_t full_step_ic0_factor_builds{};
    std::size_t half_step_ic0_factor_builds{};
    std::vector<AdaptiveHeatAttempt3D> attempts;
    // One entry per accepted macro-step; each combines the two half-step budgets.
    std::vector<TransientHeatStepDiagnostics3D> diagnostics;
    // Two successful linear solves per accepted macro-step.
    std::vector<linalg::SolverResult> substep_solvers;
};

// Adaptive Backward Euler or Crank-Nicolson using step doubling.
// The mass and stiffness operators are assembled once and reused for every
// accepted and rejected trial. The separate full-/half-step IC(0) caches
// persist across trials, with matrix-exact factor invalidation. Time-dependent
// loads and boundary data are
// evaluated at the time levels required by the selected scheme.
// The exact accepted P1 DOF vector is carried into every subsequent trial.
[[nodiscard]] AdaptiveHeat3DResult solve_adaptive_heat_3d(
    const TransientHeat3DProblem& problem,
    const AdaptiveHeat3DOptions& options = {}
);

} // namespace nova::fem
