#pragma once

#include "nova/fem/transient_heat_3d.hpp"
#include "nova/physics/electrothermal_3d.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace nova::physics {

// Experimental fixed-mesh, fixed-time-step, Backward-Euler, segregated
// electrothermal coupling. The DC electrical problem is solved at each
// thermal Picard iterate; Joule heating enters the implicit thermal step.
// Reference electrical and optional thermal material may be scaled with
// cell-average temperature using strictly positive dimensionless laws.
struct TransientElectrothermal3DProblem {
    fem::ScalarDiffusion3DProblem electrical;
    fem::TransientHeat3DProblem thermal;
    // Empty law means no temperature feedback for that material.
    std::function<double(double)> electrical_conductivity_scale;
    std::function<double(double)> thermal_conductivity_scale;
};

struct TransientElectrothermal3DOptions {
    fem::TransientHeat3DOptions thermal_time{};
    std::size_t max_coupling_iterations{40};
    double relaxation{0.7};
    double absolute_tolerance_kelvin{1.0e-8};
    double relative_tolerance{1.0e-9};
    linalg::SolverOptions electrical_solver_options{};
};

struct TransientElectrothermal3DStep {
    double time{};
    std::size_t coupling_iterations{};
    double max_fixed_point_residual_kelvin{};
    double max_normalized_residual{};
    double joule_power_watts{};
    fem::TransientHeatStepDiagnostics3D heat_budget{};
};

struct TransientElectrothermal3DResult {
    // Only fully converged steps are accepted. On nonconvergence, these
    // fields describe the last *accepted* state, not the failed trial.
    std::vector<double> temperature;
    double final_time{};
    std::vector<TransientElectrothermal3DStep> accepted_steps;
    DCJouleHeating3DResult last_electrical;
    // If false, no partially converged thermal step is committed.
    bool converged{true};
    std::size_t failed_step{};  // 1-based, 0 if all steps converged.
    double failed_max_normalized_residual{};
};

// At each timestep: solve DC conduction with sigma(T_trial), form E.J,
// solve one Backward-Euler heat step with k(T_trial), then relax the
// unrelaxed fixed-point mismatch. Dirichlet time dependence is handled by
// the underlying heat solver. No transient displacement current, Newton
// Jacobian, adaptive timestep, nonlinear capacity or radiation is modeled.
// A nonconverged coupling step returns converged=false without advancing it.
[[nodiscard]] TransientElectrothermal3DResult solve_transient_electrothermal_3d(
    const TransientElectrothermal3DProblem& problem,
    const TransientElectrothermal3DOptions& options = {}
);

} // namespace nova::physics
