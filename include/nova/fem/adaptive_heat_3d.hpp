#pragma once

#include "nova/fem/transient_heat_3d.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

// Experimental fixed-mesh adaptive stepping on the merged P1 tetrahedral
// transient solver. Compares one full step with two half steps and accepts
// the two-half-step field. All materials and boundary data remain time-independent.
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
    // Emits the Dirichlet-projected initial field (index 0) and accepted
    // macro-step fields only. Rejected/full-step trials never emit fields.
    // References are valid only for the duration of the callback.
    // Exceptions propagate; no retry occurs after an observer exception.
    std::function<void(std::size_t, double, const std::vector<double>&)> on_step;
};

struct AdaptiveHeatAttempt3D {
    double start_time{};
    double attempted_step{};
    // Maximum nodal |T_half2 - T_full| / ((2^p - 1) * (atol + rtol*scale)).
    // p=1 for backward Euler; p=2 for Crank-Nicolson. <=1 accepts.
    double normalized_error{};
    bool accepted{};
};

struct AdaptiveHeat3DResult {
    std::vector<double> temperature;
    double final_time{};
    std::size_t accepted_steps{};
    std::size_t rejected_steps{};
    std::vector<AdaptiveHeatAttempt3D> attempts;
    // One energy-budget record per accepted macro-step. Each combines
    // both half steps; heat_content is the endpoint heat content.
    std::vector<TransientHeatStepDiagnostics3D> diagnostics;
};

// Uses the current, validated fixed-step 3D heat API for every trial.
// Correctness-first port of the older adaptive #20 design: each attempt
// currently reassembles matrices (three independent solves), without
// cross-trial IC(0) or FEM-operator reuse. Not a production adaptive solver.
[[nodiscard]] AdaptiveHeat3DResult solve_adaptive_heat_3d(
    const TransientHeat3DProblem& problem,
    const AdaptiveHeat3DOptions& options = {}
);

} // namespace nova::fem
