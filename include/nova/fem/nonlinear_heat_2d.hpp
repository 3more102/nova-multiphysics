#pragma once

#include "nova/fem/transient_heat_2d.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::fem {

// Positive isotropic conductivity k(x, y, T) evaluated at the FE field.
// Capacity and boundary coefficients retain the linear thermal solver's API.
using TemperatureConductivity2D = std::function<double(double, double, double)>;

struct NonlinearHeat2DProblem {
    // Reuse source, initial state, capacity, and time-dependent edge BCs.
    // thermal.conductivity is overwritten internally on each Picard iteration.
    TransientHeat2DProblem thermal;
    TemperatureConductivity2D conductivity;
};

struct PicardHeat2DOptions {
    H1TriangleOrder order{H1TriangleOrder::P1};
    double start_time{0.0};
    double time_step{0.1};
    std::size_t steps{1U};
    std::size_t max_iterations{50U};
    // Infinity-norm update criterion:
    // max_i |u_new-u_guess| / (atol + rtol*max(|u_new|,|u_guess|)) <= 1.
    double absolute_tolerance{1.0e-10};
    double relative_tolerance{1.0e-8};
    // Under-relaxation 0 < omega <= 1; 1 = undamped Picard.
    double relaxation{1.0};
    linalg::SolverOptions solver_options{};
};

struct PicardHeatStep2D {
    double time{};
    std::size_t iterations{};
    double max_absolute_update{};
    double normalized_update{};
};

struct PicardHeat2DResult {
    // Vertex DOFs followed by P2 edge-midpoint DOFs.
    std::vector<double> temperature;
    double final_time{};
    std::vector<PicardHeatStep2D> step_diagnostics;
    // One diagnostic from the inner linear solver per Picard iteration.
    std::vector<linalg::SolverResult> linear_solvers;
};

// Fixed-time-step Backward Euler with successive conductivity freezing:
//   M (u^{n+1}-u^n)/dt + K(k(u^{n+1}))u^{n+1} = f^{n+1} + boundary terms.
//
// Each Picard trial solves the existing linear P1/P2 thermal problem, using
// k evaluated on the latest FE trial field. This is not a consistent-tangent
// Newton solver; convergence is monitored by the scaled DOF update.
// Throws on invalid inputs, nonpositive/nonfinite conductivity, inner-solver
// failure, or lack of convergence in the configured iteration budget.
// A future assembly-native implementation should avoid repeated mesh searches.
[[nodiscard]] PicardHeat2DResult solve_picard_heat_2d(
    const NonlinearHeat2DProblem& problem,
    const PicardHeat2DOptions& options = {}
);

} // namespace nova::fem
