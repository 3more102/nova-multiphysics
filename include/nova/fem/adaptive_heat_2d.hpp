#pragma once

#include "nova/fem/transient_heat_2d.hpp"

#include <cstddef>
#include <vector>

namespace nova::fem {

// Adaptive time stepping for the linear, transient P1/P2 H1 heat equation.
// A full step and two half steps share the same initial FE DOF state.
// The accepted solution is the two-half-step solution, not an extrapolation.
struct AdaptiveHeat2DOptions {
    H1TriangleOrder order{H1TriangleOrder::P1};
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

struct AdaptiveHeatAttempt2D {
    double start_time{};
    double attempted_step{};
    // RMS over all H1 DOFs: |u_2half-u_full|/(2^p-1) normalized
    // by abs_tol + rel_tol*max(|u_old|,|u_2half|).
    // A value <=1 permits acceptance.
    double normalized_error{};
    bool accepted{};
};

struct AdaptiveHeat2DResult {
    std::vector<double> temperature;
    double final_time{};
    std::size_t accepted_steps{};
    std::size_t rejected_steps{};
    std::vector<AdaptiveHeatAttempt2D> attempts;
    // Per accepted step: sums of two accepted half-step budgets.
    std::vector<TransientHeatStepDiagnostics2D> diagnostics;
    // Two successful half-step linear solver diagnostics per accepted step.
    std::vector<linalg::SolverResult> substep_solvers;
};

// Fixed-mesh adaptive time integration by step doubling.
// p=1 for Backward Euler, p=2 for Crank-Nicolson.
// Absolute tolerance must be positive; relative tolerance may be zero.
// Throws when tolerances cannot be met above min_step, when the attempt
// limit is exceeded, or when underlying numerical solves fail.
//
// Time-independent mass/stiffness operators are assembled once per adaptive
// run and reused across all trials. Time-dependent loads and dt-dependent
// systems are rebuilt; factorization/preconditioner reuse remains future work.
[[nodiscard]] AdaptiveHeat2DResult solve_adaptive_heat_2d(
    const TransientHeat2DProblem& problem,
    const AdaptiveHeat2DOptions& options = {}
);

} // namespace nova::fem
