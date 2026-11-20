#include "nova/fem/adaptive_heat_3d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {

void validate_options(const AdaptiveHeat3DOptions& options) {
    const double duration = options.end_time - options.start_time;
    if (!std::isfinite(options.start_time) || !std::isfinite(options.end_time) ||
        !std::isfinite(duration) || !(duration > 0.0)) {
        throw std::invalid_argument("Adaptive 3D heat needs finite increasing start/end times.");
    }
    if (!std::isfinite(options.initial_step) || !(options.initial_step > 0.0) ||
        !std::isfinite(options.min_step) || !(options.min_step > 0.0) ||
        !std::isfinite(options.max_step) || options.max_step < options.min_step) {
        throw std::invalid_argument("Adaptive 3D step bounds must be finite, positive and ordered.");
    }
    if (!std::isfinite(options.absolute_tolerance) || !(options.absolute_tolerance > 0.0) ||
        !std::isfinite(options.relative_tolerance) || options.relative_tolerance < 0.0) {
        throw std::invalid_argument("Adaptive 3D absolute tolerance must be positive and relative tolerance nonnegative.");
    }
    if (options.max_attempts == 0U) {
        throw std::invalid_argument("Adaptive 3D maximum attempts must be positive.");
    }
    if (options.scheme != TransientHeat3DScheme::BackwardEuler &&
        options.scheme != TransientHeat3DScheme::CrankNicolson) {
        throw std::invalid_argument("Adaptive 3D time scheme is unsupported.");
    }
}

[[nodiscard]] double normalized_error(
    const std::vector<double>& previous,
    const std::vector<double>& full,
    const std::vector<double>& half,
    const AdaptiveHeat3DOptions& options
) {
    if (previous.empty() || previous.size() != full.size() || full.size() != half.size()) {
        throw std::invalid_argument("Adaptive 3D step trial field sizes differ.");
    }
    const double divisor = options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 3.0;
    double norm = 0.0;
    for (std::size_t i = 0; i < previous.size(); ++i) {
        const double scale = options.absolute_tolerance +
            options.relative_tolerance * std::max(std::abs(previous[i]), std::abs(half[i]));
        if (!std::isfinite(scale) || !(scale > 0.0)) {
            throw std::overflow_error("Adaptive 3D error scale is non-finite.");
        }
        const double component = ((half[i] - full[i]) / divisor) / scale;
        if (!std::isfinite(component)) {
            throw std::overflow_error("Adaptive 3D error component is non-finite.");
        }
        // hypot avoids overflow/underflow of naive sum(component*component).
        norm = std::hypot(norm, component);
    }
    return norm / std::sqrt(static_cast<double>(previous.size()));
}

[[nodiscard]] TransientHeatStepDiagnostics3D combined_budget(
    const std::vector<TransientHeatStepDiagnostics3D>& half_steps
) {
    if (half_steps.size() != 2U) {
        throw std::logic_error("Adaptive 3D trial needs exactly two half-step budgets.");
    }
    auto budget = half_steps.back();
    budget.source_energy += half_steps.front().source_energy;
    budget.neumann_energy += half_steps.front().neumann_energy;
    budget.robin_energy += half_steps.front().robin_energy;
    budget.dirichlet_energy += half_steps.front().dirichlet_energy;
    budget.balance_residual += half_steps.front().balance_residual;
    return budget;
}

} // namespace

AdaptiveHeat3DResult solve_adaptive_heat_3d(
    const TransientHeat3DProblem& problem,
    const AdaptiveHeat3DOptions& options
) {
    validate_options(options);

    // Keep the exact accepted nodal state as the restart for both trial paths.
    // Existing transient solve validates and projects prescribed nodes.
    TransientHeat3DProblem trial = problem;
    std::vector<double> state = problem.initial_temperature;
    AdaptiveHeat3DResult result;
    double time = options.start_time;
    double step = std::clamp(options.initial_step, options.min_step, options.max_step);
    const double order = options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 2.0;

    while (time < options.end_time) {
        if (result.attempts.size() >= options.max_attempts) {
            throw std::runtime_error("Adaptive 3D heat maximum attempts exhausted.");
        }
        const double remaining = options.end_time - time;
        const double dt = std::min(step, remaining);
        if (!std::isfinite(dt) || !(dt > 0.0) || !(time + dt > time) ||
            !(time + 0.5 * dt > time)) {
            throw std::runtime_error("Adaptive 3D time step cannot advance representable time.");
        }

        trial.initial_temperature = state;
        TransientHeat3DOptions fixed;
        fixed.scheme = options.scheme;
        fixed.start_time = time;
        fixed.time_step = dt;
        fixed.steps = 1U;
        fixed.solver_options = options.solver_options;
        fixed.retain_step_solutions = false;
        const auto full = solve_transient_heat_3d(trial, fixed);

        fixed.time_step = dt * 0.5;
        fixed.steps = 2U;
        auto half = solve_transient_heat_3d(trial, fixed);

        const double error = normalized_error(state, full.temperature, half.temperature, options);
        const bool accepted = error <= 1.0;
        result.attempts.push_back({time, dt, error, accepted});

        if (accepted) {
            state = std::move(half.temperature);
            result.diagnostics.push_back(combined_budget(half.diagnostics));
            for (auto& solver : half.linear_solvers) {
                result.substep_solvers.push_back(std::move(solver));
            }
            ++result.accepted_steps;
            time = std::min(options.end_time, time + dt);
        } else {
            ++result.rejected_steps;
            if (dt <= options.min_step) {
                throw std::runtime_error("Adaptive 3D error exceeds tolerance at minimum step.");
            }
        }

        const double bounded_error = std::max(error, 1.0e-14);
        const double factor = std::clamp(
            0.9 * std::pow(bounded_error, -1.0 / (order + 1.0)),
            0.2, accepted ? 4.0 : 0.9);
        step = std::clamp(dt * factor, options.min_step, options.max_step);
    }

    result.temperature = std::move(state);
    result.final_time = time;
    return result;
}

} // namespace nova::fem
