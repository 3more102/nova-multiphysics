#include "nova/fem/adaptive_heat_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

void validate_options(const AdaptiveHeat3DOptions& options) {
    if (!std::isfinite(options.start_time) ||
        !std::isfinite(options.end_time) ||
        !(options.end_time > options.start_time)) {
        throw std::invalid_argument("Adaptive 3D heat requires finite increasing times.");
    }
    if (!std::isfinite(options.initial_step) || !(options.initial_step > 0.0) ||
        !std::isfinite(options.min_step) || !(options.min_step > 0.0) ||
        !std::isfinite(options.max_step) || options.max_step < options.min_step) {
        throw std::invalid_argument("Adaptive 3D time-step bounds are invalid.");
    }
    if (!std::isfinite(options.absolute_tolerance) ||
        !(options.absolute_tolerance > 0.0) ||
        !std::isfinite(options.relative_tolerance) ||
        options.relative_tolerance < 0.0) {
        throw std::invalid_argument("Adaptive 3D tolerances are invalid.");
    }
    if (options.max_attempts == 0U) {
        throw std::invalid_argument("Adaptive 3D max_attempts must be positive.");
    }
    if (options.scheme != TransientHeat3DScheme::BackwardEuler &&
        options.scheme != TransientHeat3DScheme::CrankNicolson) {
        throw std::invalid_argument("Adaptive 3D time scheme is unsupported.");
    }
}

double normalized_step_error(
    const std::vector<double>& initial,
    const std::vector<double>& full,
    const std::vector<double>& half,
    const AdaptiveHeat3DOptions& options
) {
    if (initial.empty() || full.size() != initial.size() ||
        half.size() != initial.size()) {
        throw std::invalid_argument("Adaptive 3D trial state dimensions disagree.");
    }

    const double divisor =
        options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 3.0;
    long double squared = 0.0L;
    for (std::size_t i = 0; i < initial.size(); ++i) {
        const double scale = options.absolute_tolerance +
            options.relative_tolerance *
            std::max(std::abs(initial[i]), std::abs(half[i]));
        if (!(scale > 0.0) || !std::isfinite(scale)) {
            throw std::overflow_error("Adaptive 3D error scale is non-finite.");
        }
        const double difference = half[i] - full[i];
        const double error = (difference / divisor) / scale;
        if (!std::isfinite(error)) {
            throw std::overflow_error("Adaptive 3D local error is non-finite.");
        }
        squared += static_cast<long double>(error) * error;
    }
    const long double rms = std::sqrt(
        squared / static_cast<long double>(initial.size()));
    if (!std::isfinite(rms) ||
        rms > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Adaptive 3D RMS error overflowed.");
    }
    return static_cast<double>(rms);
}

TransientHeatStepDiagnostics3D combine_half_budgets(
    const std::vector<TransientHeatStepDiagnostics3D>& half
) {
    if (half.size() != 2U) {
        throw std::logic_error("Adaptive 3D expected two half-step budgets.");
    }
    TransientHeatStepDiagnostics3D combined = half.back();
    combined.source_energy += half.front().source_energy;
    combined.neumann_energy += half.front().neumann_energy;
    combined.robin_energy += half.front().robin_energy;
    combined.dirichlet_energy += half.front().dirichlet_energy;
    combined.balance_residual += half.front().balance_residual;
    return combined;
}

} // namespace

AdaptiveHeat3DResult solve_adaptive_heat_3d(
    const TransientHeat3DProblem& problem,
    const AdaptiveHeat3DOptions& options
) {
    validate_options(options);
    if (problem.initial_temperature.size() != problem.diffusion.mesh.nodes.size()) {
        throw std::invalid_argument("Adaptive 3D initial field length mismatch.");
    }
    for (const double value : problem.initial_temperature) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Adaptive 3D initial temperature must be finite.");
        }
    }

    auto trial_problem = problem;
    std::vector<double> state = problem.initial_temperature;
    AdaptiveHeat3DResult result;
    double time = options.start_time;
    double step = std::clamp(
        options.initial_step, options.min_step, options.max_step);
    const double order =
        options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 2.0;
    bool initial_state_emitted = false;

    while (time < options.end_time) {
        if (result.attempts.size() >= options.max_attempts) {
            throw std::runtime_error("Adaptive 3D maximum attempts exceeded.");
        }
        const double dt = std::min(step, options.end_time - time);
        if (!(dt > 0.0) || !std::isfinite(dt) ||
            !std::isfinite(time + dt) || !(time + dt > time)) {
            throw std::runtime_error("Adaptive 3D time cannot advance representably.");
        }

        trial_problem.initial_temperature = state;
        TransientHeat3DOptions fixed;
        fixed.scheme = options.scheme;
        fixed.start_time = time;
        fixed.time_step = dt;
        fixed.steps = 1U;
        fixed.solver_options = options.solver_options;
        fixed.retain_step_solutions = false;
        // The first full-step solve supplies the validated Dirichlet projection;
        // its coarse step result is deliberately not sent to the user observer.
        if (!initial_state_emitted) {
            fixed.on_step = [&](std::size_t step_index, double physical_time,
                                const std::vector<double>& projected_state) {
                if (step_index == 0U) {
                    state = projected_state;
                    initial_state_emitted = true;
                    if (options.on_step) options.on_step(0U, physical_time, state);
                }
            };
        }
        const auto full = solve_transient_heat_3d(trial_problem, fixed);
        fixed.on_step = {};
        // Carry the projected initial state into both half-step trials.
        trial_problem.initial_temperature = state;

        fixed.time_step = 0.5 * dt;
        fixed.steps = 2U;
        auto half = solve_transient_heat_3d(trial_problem, fixed);
        const double error = normalized_step_error(
            state, full.temperature, half.temperature, options);
        const bool accepted = error <= 1.0;
        result.attempts.push_back({time, dt, error, accepted});

        if (accepted) {
            state = std::move(half.temperature);
            result.diagnostics.push_back(combine_half_budgets(half.diagnostics));
            for (auto& solver : half.linear_solvers) {
                result.substep_solvers.push_back(std::move(solver));
            }
            ++result.accepted_steps;
            time = std::min(options.end_time, time + dt);
            if (options.on_step) {
                options.on_step(result.accepted_steps, time, state);
            }
        } else {
            ++result.rejected_steps;
            if (dt <= options.min_step) {
                throw std::runtime_error(
                    "Adaptive 3D error exceeds tolerance at minimum time step.");
            }
        }

        const double bounded_error = std::max(error, 1.0e-14);
        const double proposed_factor =
            0.9 * std::pow(bounded_error, -1.0 / (order + 1.0));
        const double factor =
            std::clamp(proposed_factor, 0.2, accepted ? 4.0 : 0.9);
        step = std::clamp(dt * factor, options.min_step, options.max_step);
    }

    result.temperature = std::move(state);
    result.final_time = time;
    return result;
}

} // namespace nova::fem
