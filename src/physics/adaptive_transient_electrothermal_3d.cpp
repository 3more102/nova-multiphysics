#include "nova/physics/adaptive_transient_electrothermal_3d.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>

namespace nova::physics {
namespace {

void validate_progress(const AdaptiveTransientElectrothermal3DProgress& progress) {
    if (progress.accepted_steps > progress.attempted_steps ||
        progress.rejected_steps !=
            progress.attempted_steps - progress.accepted_steps) {
        throw std::invalid_argument("Invalid cumulative electrothermal attempt counts.");
    }
    for (const double value : {progress.source_energy, progress.neumann_energy,
                               progress.robin_energy, progress.dirichlet_energy,
                               progress.balance_residual, progress.latest_heat_content}) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Nonfinite cumulative electrothermal energy.");
        }
    }
}

void advance_attempt(AdaptiveTransientElectrothermal3DProgress& progress) {
    if (progress.attempted_steps == std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error("Cumulative electrothermal attempt count overflows.");
    }
    ++progress.attempted_steps;
}

void accumulate_accepted_step(
    AdaptiveTransientElectrothermal3DProgress& progress,
    const AdaptiveTransientElectrothermal3DStep& step
) {
    if (progress.accepted_steps == std::numeric_limits<std::size_t>::max() ||
        progress.nonlinear_iterations >
            std::numeric_limits<std::size_t>::max() - step.nonlinear_iterations) {
        throw std::overflow_error("Cumulative electrothermal counters overflow.");
    }
    auto next = progress;
    ++next.accepted_steps;
    next.nonlinear_iterations += step.nonlinear_iterations;
    const auto add = [](double lhs, double rhs) {
        const double sum = lhs + rhs;
        if (!std::isfinite(sum)) {
            throw std::overflow_error("Cumulative electrothermal energy overflows.");
        }
        return sum;
    };
    next.source_energy = add(next.source_energy, step.heat_budget.source_energy);
    next.neumann_energy = add(next.neumann_energy, step.heat_budget.neumann_energy);
    next.robin_energy = add(next.robin_energy, step.heat_budget.robin_energy);
    next.dirichlet_energy = add(next.dirichlet_energy, step.heat_budget.dirichlet_energy);
    next.balance_residual = add(next.balance_residual,
                                step.heat_budget.balance_residual);
    next.latest_heat_content = step.heat_budget.heat_content;
    progress = next;
}


void validate(const TransientElectrothermal3DProblem& problem,
              const AdaptiveTransientElectrothermal3DOptions& o) {
    const double start = o.coupling.thermal_time.start_time;
    const auto positive_finite = [](double x) {
        return std::isfinite(x) && x > 0.0;
    };
    if (!std::isfinite(start) || !std::isfinite(o.end_time) ||
        !(o.end_time > start) ||
        !positive_finite(o.initial_time_step) ||
        !positive_finite(o.minimum_time_step) ||
        !positive_finite(o.maximum_time_step) ||
        !(o.minimum_time_step <= o.initial_time_step) ||
        !(o.initial_time_step <= o.maximum_time_step) ||
        !std::isfinite(o.absolute_error_kelvin) ||
        o.absolute_error_kelvin < 0.0 ||
        !std::isfinite(o.relative_error) || o.relative_error < 0.0 ||
        (o.absolute_error_kelvin == 0.0 && o.relative_error == 0.0) ||
        !positive_finite(o.safety_factor) || o.safety_factor > 1.0 ||
        !positive_finite(o.minimum_shrink_factor) ||
        o.minimum_shrink_factor >= 1.0 ||
        !std::isfinite(o.maximum_growth_factor) ||
        o.maximum_growth_factor <= 1.0 ||
        o.maximum_attempts == 0U ||
        o.maximum_accepted_steps == 0U ||
        o.coupling.thermal_time.scheme != fem::TransientHeat3DScheme::BackwardEuler ||
        o.coupling.thermal_time.on_step) {
        throw std::invalid_argument("Invalid adaptive electrothermal step-control options.");
    }
    validate_progress(o.previous_progress);
    if (problem.thermal.initial_temperature.size() !=
        problem.electrical.mesh.nodes.size()) {
        throw std::invalid_argument("Adaptive electrothermal initial field size mismatch.");
    }
}

fem::TransientHeatStepDiagnostics3D sum_half_step_budgets(
    const fem::TransientHeatStepDiagnostics3D& first,
    const fem::TransientHeatStepDiagnostics3D& second
) {
    fem::TransientHeatStepDiagnostics3D combined;
    combined.time = second.time;
    combined.heat_content = second.heat_content;
    combined.source_energy = first.source_energy + second.source_energy;
    combined.neumann_energy = first.neumann_energy + second.neumann_energy;
    combined.robin_energy = first.robin_energy + second.robin_energy;
    combined.dirichlet_energy = first.dirichlet_energy + second.dirichlet_energy;
    combined.balance_residual = first.balance_residual + second.balance_residual;
    for (const auto value : {combined.heat_content, combined.source_energy,
                             combined.neumann_energy, combined.robin_energy,
                             combined.dirichlet_energy, combined.balance_residual}) {
        if (!std::isfinite(value)) {
            throw std::overflow_error("Adaptive electrothermal accumulated energy overflows.");
        }
    }
    return combined;
}

TransientElectrothermal3DResult solve_one_step(
    const TransientElectrothermal3DProblem& problem,
    const TransientElectrothermal3DOptions& base,
    const std::vector<double>& accepted_temperature,
    double start_time,
    double step_size
) {
    auto p = problem;
    p.thermal.initial_temperature = accepted_temperature;
    auto o = base;
    o.thermal_time.start_time = start_time;
    o.thermal_time.time_step = step_size;
    o.thermal_time.steps = 1U;
    o.thermal_time.retain_step_solutions = false;
    return solve_transient_electrothermal_3d(p, o);
}

double estimate_normalized_error(
    const std::vector<double>& one_full,
    const std::vector<double>& two_half,
    const AdaptiveTransientElectrothermal3DOptions& o
) {
    if (one_full.size() != two_half.size() || one_full.empty()) {
        throw std::runtime_error("Adaptive electrothermal candidate field size mismatch.");
    }
    double error = 0.0;
    for (std::size_t i = 0; i < one_full.size(); ++i) {
        const double a = one_full[i], b = two_half[i];
        if (!std::isfinite(a) || !std::isfinite(b)) {
            throw std::runtime_error("Adaptive electrothermal trial temperature is nonfinite.");
        }
        const double scale = o.absolute_error_kelvin +
            o.relative_error * std::max(std::abs(a), std::abs(b));
        const double e = std::abs(a - b) / scale;
        if (!std::isfinite(e)) {
            throw std::overflow_error("Adaptive electrothermal local error overflows.");
        }
        error = std::max(error, e);
    }
    return error;
}

double next_time_step(double current, double error,
                      const AdaptiveTransientElectrothermal3DOptions& o,
                      bool accepted, bool nonlinear_converged) {
    double factor = o.minimum_shrink_factor;
    if (!nonlinear_converged) {
        factor = std::max(o.minimum_shrink_factor, 0.5);
    } else if (error == 0.0) {
        factor = o.maximum_growth_factor;
    } else {
        factor = o.safety_factor / std::sqrt(error);
        factor = std::clamp(factor, o.minimum_shrink_factor, o.maximum_growth_factor);
    }
    if (!accepted) factor = std::min(factor, 0.9);
    const double scaled = static_cast<double>(
        static_cast<long double>(current) * factor);
    if (!std::isfinite(scaled)) return o.maximum_time_step;
    return std::clamp(scaled, o.minimum_time_step, o.maximum_time_step);
}

} // namespace

AdaptiveTransientElectrothermal3DResult solve_adaptive_transient_electrothermal_3d(
    const TransientElectrothermal3DProblem& problem,
    const AdaptiveTransientElectrothermal3DOptions& options
) {
    validate(problem, options);
    AdaptiveTransientElectrothermal3DResult result;
    result.progress = options.previous_progress;
    result.temperature = problem.thermal.initial_temperature;
    result.final_time = options.coupling.thermal_time.start_time;
    result.attempts.reserve(std::min<std::size_t>(options.maximum_attempts, 10000U));
    result.accepted_steps.reserve(
        std::min<std::size_t>(options.maximum_accepted_steps, 5000U));
    double dt = options.initial_time_step;

    while (result.final_time < options.end_time) {
        if (result.attempts.size() >= options.maximum_attempts) {
            result.stop_reason = AdaptiveTransientElectrothermal3DStopReason::AttemptLimit;
            return result;
        }
        if (result.accepted_steps.size() >= options.maximum_accepted_steps) {
            result.stop_reason = AdaptiveTransientElectrothermal3DStopReason::AcceptedStepLimit;
            return result;
        }

        const double start = result.final_time;
        const double remaining = options.end_time - start;
        const double h = std::min(dt, remaining);
        const double end = start + h;
        const double half = h / 2.0;
        const double middle = start + half;
        if (!std::isfinite(h) || h <= 0.0 ||
            !(middle > start) || !(end > middle) ||
            end > options.end_time) {
            result.stop_reason = AdaptiveTransientElectrothermal3DStopReason::TimeResolution;
            return result;
        }
        advance_attempt(result.progress);
        AdaptiveTransientElectrothermal3DAttempt attempt;
        attempt.start_time = start;
        attempt.attempted_step = h;
        attempt.end_time = end;
        attempt.max_normalized_local_error = std::numeric_limits<double>::infinity();

        // All three candidates originate from the same LAST accepted state.
        // Neither the full step nor any failed half step can commit its field.
        auto full = solve_one_step(
            problem, options.coupling, result.temperature, start, h);
        if (!full.converged) {
            result.attempts.push_back(attempt);
        } else {
            auto first = solve_one_step(
                problem, options.coupling, result.temperature, start, half);
            if (!first.converged) {
                result.attempts.push_back(attempt);
            } else {
                auto second = solve_one_step(
                    problem, options.coupling, first.temperature, middle, half);
                if (!second.converged) {
                    result.attempts.push_back(attempt);
                } else {
                    attempt.nonlinear_converged = true;
                    attempt.max_normalized_local_error =
                        estimate_normalized_error(
                            full.temperature, second.temperature, options);
                    attempt.accepted = attempt.max_normalized_local_error <= 1.0;
                    result.attempts.push_back(attempt);
                    if (attempt.accepted) {
                        if (first.accepted_steps.size() != 1U ||
                            second.accepted_steps.size() != 1U ||
                            full.accepted_steps.size() != 1U) {
                            throw std::runtime_error(
                                "Adaptive electrothermal step result is inconsistent.");
                        }
                        const auto& d1 = first.accepted_steps.front();
                        const auto& d2 = second.accepted_steps.front();
                        AdaptiveTransientElectrothermal3DStep accepted_step;
                        accepted_step.time = end;
                        accepted_step.time_step = h;
                        accepted_step.max_normalized_local_error =
                            attempt.max_normalized_local_error;
                        accepted_step.next_time_step = next_time_step(
                            h, attempt.max_normalized_local_error, options, true, true);
                        accepted_step.nonlinear_iterations =
                            full.accepted_steps.front().coupling_iterations +
                            d1.coupling_iterations + d2.coupling_iterations;
                        accepted_step.joule_power_watts = d2.joule_power_watts;
                        accepted_step.heat_budget =
                            sum_half_step_budgets(d1.heat_budget, d2.heat_budget);

                        accumulate_accepted_step(result.progress, accepted_step);
                        accepted_step.progress = result.progress;
                        result.temperature = std::move(second.temperature);
                        result.last_electrical = std::move(second.last_electrical);
                        result.final_time = end;
                        result.accepted_steps.push_back(accepted_step);
                        dt = accepted_step.next_time_step;
                        // The only side-effecting observer is reached after
                        // acceptance; discarded nonlinear candidates cannot
                        // trigger a checkpoint or external output.
                        if (options.on_accepted_step) {
                            options.on_accepted_step(
                                result.accepted_steps.back(), result.temperature);
                        }
                        continue;
                    }
                }
            }
        }
        // Failed nonlinear or temporal error-control trial: roll back.
        ++result.rejected_steps;
        if (result.progress.rejected_steps ==
            std::numeric_limits<std::size_t>::max()) {
            throw std::overflow_error("Cumulative electrothermal rejection count overflows.");
        }
        ++result.progress.rejected_steps;
        if (h <= options.minimum_time_step) {
            result.stop_reason =
                AdaptiveTransientElectrothermal3DStopReason::MinimumTimeStep;
            return result;
        }
        const auto& last = result.attempts.back();
        dt = next_time_step(h, last.max_normalized_local_error, options,
                            false, last.nonlinear_converged);
    }
    result.completed = true;
    result.stop_reason = AdaptiveTransientElectrothermal3DStopReason::Completed;
    return result;
}

} // namespace nova::physics
