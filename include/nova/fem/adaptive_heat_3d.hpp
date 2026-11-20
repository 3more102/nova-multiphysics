#pragma once

#include "nova/fem/transient_heat_3d.hpp"

#include <cstddef>
#include <vector>

namespace nova::fem {

// Fixed-mesh adaptive stepping for the existing time-independent-coefficient
// P1 tetrahedral transient heat solver. Compare one full step with two half
// steps; only the two-half-step solution is accepted.
struct AdaptiveHeat3DOptions {
    TransientHeat3DScheme scheme{TransientHeat3DScheme::CrankNicolson};
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
    // RMS of |T_two_half-T_full|/(2^p-1), normalized by
    // atol+rtol*max(|T_old|,|T_two_half|), p=1 (BE) or 2 (CN).
    double normalized_error{};
    bool accepted{};
};

struct AdaptiveHeat3DResult {
    std::vector<double> temperature;
    double final_time{};
    std::size_t accepted_steps{};
    std::size_t rejected_steps{};
    std::vector<AdaptiveHeatAttempt3D> attempts;
    // One heat budget for each accepted step, summing both half steps.
    std::vector<TransientHeatStepDiagnostics3D> diagnostics;
    // Two solver diagnostic records per accepted step; solutions are omitted
    // to avoid duplicating large field vectors in the result.
    std::vector<linalg::SolverResult> substep_solvers;
};

// Mainline integration of the adaptive controller from draft PR #20.
// Uses the current validated fixed-step solver and its numerical diagnostics.
// Each trial currently assembles FEM operators anew; shared prepared operators
// and cross-trial IC(0) cache reuse are separate follow-up optimizations.
// Throws for invalid controls, solver failure, attempts exhaustion, or an
// unachievable error tolerance at the minimum step.
[[nodiscard]] AdaptiveHeat3DResult solve_adaptive_heat_3d(
    const TransientHeat3DProblem& problem,
    const AdaptiveHeat3DOptions& options = {}
);

} // namespace nova::fem
