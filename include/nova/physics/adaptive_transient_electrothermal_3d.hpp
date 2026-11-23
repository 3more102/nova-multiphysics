#pragma once

#include "nova/physics/transient_electrothermal_3d.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::physics {

struct AdaptiveTransientElectrothermal3DStep;

// Experimental adaptive-step, fixed-mesh, Backward Euler electrothermal
// coupling. Each trial compares one full step to two independently
// converged half steps. Accepted state uses the TWO half steps.
struct AdaptiveTransientElectrothermal3DOptions {
    TransientElectrothermal3DOptions coupling{};
    double end_time{1.0};
    double initial_time_step{0.1};
    double minimum_time_step{1.0e-6};
    double maximum_time_step{0.5};
    double absolute_error_kelvin{1.0e-4};
    double relative_error{1.0e-6};
    double safety_factor{0.9};          // (0,1]
    double minimum_shrink_factor{0.2};  // (0,1)
    double maximum_growth_factor{2.0}; // >1
    std::size_t maximum_attempts{10000U};
    std::size_t maximum_accepted_steps{5000U};
    // Invoked exactly once AFTER an accepted step is committed, never for
    // full-step comparisons or rejected trials. The temperature is the
    // accepted two-half-step field, and can be serialized with the existing
    // mesh-validated transient heat checkpoint writer. Exceptions propagate.
    std::function<void(const AdaptiveTransientElectrothermal3DStep&,
                       const std::vector<double>&)> on_accepted_step;
};

enum class AdaptiveTransientElectrothermal3DStopReason {
    Completed,
    AttemptLimit,
    AcceptedStepLimit,
    MinimumTimeStep,
    TimeResolution
};

struct AdaptiveTransientElectrothermal3DAttempt {
    double start_time{};
    double attempted_step{};
    double end_time{};
    // infinity denotes that one of three nonlinear solves did not converge.
    double max_normalized_local_error{};
    bool nonlinear_converged{};
    bool accepted{};
};

struct AdaptiveTransientElectrothermal3DStep {
    double time{};
    double time_step{};
    double max_normalized_local_error{};
    // Controller's proposed dt for the NEXT attempt, clamped to the
    // configured [minimum_time_step, maximum_time_step]. Save alongside
    // the accepted time/temperature for deterministic step-boundary restart.
    double next_time_step{};
    // Sum of all nonlinear iterations over the full/half-step comparison,
    // including the discarded full-step candidate.
    std::size_t nonlinear_iterations{};
    double joule_power_watts{}; // power at the final accepted half-step
    // Conservative sum of two accepted half-step discrete energy budgets.
    // heat_content is from the second half-step; all energy inputs and
    // balance residuals are summed, and time is the accepted end time.
    fem::TransientHeatStepDiagnostics3D heat_budget{};
};

// On early termination, temperature and final_time refer ONLY to the
// last accepted end state. Exceptions from invalid physics and failed
// inner linear systems propagate; nonlinear nonconvergence is rejected.
struct AdaptiveTransientElectrothermal3DResult {
    std::vector<double> temperature;
    double final_time{};
    bool completed{};
    AdaptiveTransientElectrothermal3DStopReason stop_reason{
        AdaptiveTransientElectrothermal3DStopReason::AttemptLimit};
    std::vector<AdaptiveTransientElectrothermal3DAttempt> attempts;
    std::vector<AdaptiveTransientElectrothermal3DStep> accepted_steps;
    DCJouleHeating3DResult last_electrical;
    std::size_t rejected_steps{};
};

// The local temporal error indicator is the maximal normalized nodal
// difference between the one-full and two-half Backward Euler solutions,
// with absolute+relative Kelvin scaling. This is not a certified continuum
// error estimate or an automatic spatial mesh refinement method.
[[nodiscard]] AdaptiveTransientElectrothermal3DResult
solve_adaptive_transient_electrothermal_3d(
    const TransientElectrothermal3DProblem& problem,
    const AdaptiveTransientElectrothermal3DOptions& options = {}
);

} // namespace nova::physics
