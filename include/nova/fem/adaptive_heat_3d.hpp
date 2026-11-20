#pragma once

#include "nova/fem/transient_heat_3d.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

// Experimental fixed-mesh temporal adaptivity for the existing P1 tetrahedral
// solver. All coefficients and boundary data remain time-independent.
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
    // macro-step fields only. Rejected/full-step trials are never published.
    // The vector reference is valid only during the callback; exceptions propagate.
    std::function<void(std::size_t, double, const std::vector<double>&)> on_step;
};

struct AdaptiveHeatAttempt3D {
    double start_time{};
    double attempted_step{};
    // RMS over all nodal DOFs of the Richardson-scaled local error,
    // normalized with abs_tol + rel_tol*max(|T_old|, |T_2half|).
    double normalized_error{};
    bool accepted{};
};

struct AdaptiveHeat3DResult {
    std::vector<double> temperature;
    double final_time{};
    std::size_t accepted_steps{};
    std::size_t rejected_steps{};
    std::vector<AdaptiveHeatAttempt3D> attempts;
    // Each entry combines the two accepted half-step energy budgets.
    std::vector<TransientHeatStepDiagnostics3D> diagnostics;
    // Two solver diagnostics per accepted time step; temperature vectors
    // are omitted because the full adaptive field history is not retained.
    std::vector<linalg::SolverResult> substep_solvers;
};

// Compare a full step with two half steps from the same initial nodal state.
// Accept the two-half-step field, not an extrapolation. Rejected attempts
// never modify the accepted field or its energy diagnostics.
// p=1 for backward Euler, p=2 for Crank-Nicolson.
// Time-independent FEM operators are assembled once per integration. Two
// bounded dt-specific caches retain the full- and half-step systems; a changed
// step size or scheme rebuilds the affected system. IC(0) factors are reused
// only while the exact constrained CSR matrix remains unchanged.
// Throws on invalid options, solver failure, unreachable tolerance at
// min_step, non-advancing floating-point time, or max_attempts exhaustion.
[[nodiscard]] AdaptiveHeat3DResult solve_adaptive_heat_3d(
    const TransientHeat3DProblem& problem,
    const AdaptiveHeat3DOptions& options = {}
);

} // namespace nova::fem
