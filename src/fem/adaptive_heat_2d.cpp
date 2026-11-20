#include "nova/fem/adaptive_heat_2d.hpp"

#include "transient_heat_2d_internal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

void validate_options(const AdaptiveHeat2DOptions& options) {
    if (!std::isfinite(options.start_time) || !std::isfinite(options.end_time) ||
        options.end_time <= options.start_time) {
        throw std::invalid_argument("Adaptive heat needs finite increasing start/end times.");
    }
    if (!std::isfinite(options.initial_step) || options.initial_step <= 0.0 ||
        !std::isfinite(options.min_step) || options.min_step <= 0.0 ||
        !std::isfinite(options.max_step) || options.max_step < options.min_step) {
        throw std::invalid_argument("Adaptive time-step bounds must be finite, positive and ordered.");
    }
    if (!std::isfinite(options.absolute_tolerance) ||
        options.absolute_tolerance <= 0.0 ||
        !std::isfinite(options.relative_tolerance) ||
        options.relative_tolerance < 0.0) {
        throw std::invalid_argument("Adaptive absolute tolerance must be positive; relative tolerance nonnegative.");
    }
    if (options.max_attempts == 0U) {
        throw std::invalid_argument("Adaptive maximum attempts must be positive.");
    }
    if (options.scheme != TransientTimeScheme::BackwardEuler &&
        options.scheme != TransientTimeScheme::CrankNicolson) {
        throw std::invalid_argument("Adaptive thermal time scheme is unsupported.");
    }
}

[[nodiscard]] double normalized_step_error(
    const std::vector<double>& initial,
    const std::vector<double>& full,
    const std::vector<double>& half,
    const AdaptiveHeat2DOptions& options
) {
    if (full.size() != half.size() || full.size() != initial.size() ||
        initial.empty()) {
        throw std::invalid_argument("Adaptive trial field dimension mismatch.");
    }
    const double divisor =
        options.scheme == TransientTimeScheme::BackwardEuler ? 1.0 : 3.0;
    double sum = 0.0;
    for (std::size_t i = 0; i < full.size(); ++i) {
        const double scale = options.absolute_tolerance +
            options.relative_tolerance * std::max(
                std::abs(initial[i]), std::abs(half[i])
            );
        const double error = ((half[i] - full[i]) / divisor) / scale;
        if (!std::isfinite(error)) {
            throw std::overflow_error("Adaptive local truncation error is non-finite.");
        }
        sum += error * error;
        if (!std::isfinite(sum)) {
            throw std::overflow_error("Adaptive squared error accumulation overflowed.");
        }
    }
    return std::sqrt(sum / static_cast<double>(full.size()));
}

[[nodiscard]] TransientHeatStepDiagnostics2D combine_half_budgets(
    const std::vector<TransientHeatStepDiagnostics2D>& half
) {
    if (half.size() != 2U) {
        throw std::logic_error("Adaptive heat expected exactly two half-step budgets.");
    }
    TransientHeatStepDiagnostics2D diagnostic = half.back();
    diagnostic.source_energy += half.front().source_energy;
    diagnostic.neumann_energy += half.front().neumann_energy;
    diagnostic.robin_energy += half.front().robin_energy;
    diagnostic.dirichlet_energy += half.front().dirichlet_energy;
    diagnostic.balance_residual += half.front().balance_residual;
    return diagnostic;
}

} // namespace

AdaptiveHeat2DResult solve_adaptive_heat_2d(
    const TransientHeat2DProblem& problem,
    const AdaptiveHeat2DOptions& options
) {
    validate_options(options);

    const H1FunctionSpace2D space(problem.mesh, options.order);
    // Material matrices, DOF boundary topology and capacity weights are
    // prepared once. Time-dependent source and boundary data remain evaluated
    // at each required physical time level by the prepared step solver.
    std::vector<double> state = detail::initial_heat_state(problem, space);
    const auto prepared = detail::prepare_transient_heat_2d(problem, space);
    detail::TransientHeat2DStepCache full_step_cache;
    detail::TransientHeat2DStepCache half_step_cache;

    AdaptiveHeat2DResult result;
    double time = options.start_time;
    double step = std::clamp(
        options.initial_step, options.min_step, options.max_step
    );
    const double order =
        options.scheme == TransientTimeScheme::BackwardEuler ? 1.0 : 2.0;

    while (time < options.end_time) {
        if (result.attempts.size() >= options.max_attempts) {
            throw std::runtime_error("Adaptive heat maximum step attempts exhausted.");
        }
        const double remaining = options.end_time - time;
        const double dt = std::min(step, remaining);
        if (dt <= 0.0 || time + dt <= time) {
            throw std::runtime_error("Adaptive time step is too small to advance floating-point time.");
        }

        TransientHeat2DOptions fixed;
        fixed.order = options.order;
        fixed.scheme = options.scheme;
        fixed.start_time = time;
        fixed.time_step = dt;
        fixed.steps = 1U;
        fixed.solver_options = options.solver_options;
        const TransientHeat2DResult full =
            detail::solve_prepared_transient_heat_2d(
                problem, space, prepared, fixed, state, full_step_cache
            );

        fixed.time_step = 0.5 * dt;
        fixed.steps = 2U;
        TransientHeat2DResult half =
            detail::solve_prepared_transient_heat_2d(
                problem, space, prepared, fixed, state, half_step_cache
            );

        const double error = normalized_step_error(
            state, full.temperature, half.temperature, options
        );
        const bool accepted = error <= 1.0;
        result.attempts.push_back({time, dt, error, accepted});

        if (accepted) {
            // Preserve accepted DOFs and two substep diagnostics from the
            // accurate two-half-step trial. Rejected trials have no side effects.
            state = std::move(half.temperature);
            result.diagnostics.push_back(
                combine_half_budgets(half.diagnostics)
            );
            for (auto& solver : half.linear_solvers) {
                result.substep_solvers.push_back(std::move(solver));
            }
            ++result.accepted_steps;
            time = std::min(options.end_time, time + dt);
        } else {
            ++result.rejected_steps;
            if (dt <= options.min_step) {
                throw std::runtime_error(
                    "Adaptive thermal error exceeds tolerance at minimum time step."
                );
            }
        }

        const double bounded_error = std::max(error, 1.0e-14);
        double factor = 0.9 * std::pow(bounded_error, -1.0 / (order + 1.0));
        factor = std::clamp(factor, 0.2, accepted ? 4.0 : 0.9);
        step = std::clamp(dt * factor, options.min_step, options.max_step);
    }

    result.full_step_base_operator_builds = full_step_cache.base_operator_builds;
    result.half_step_base_operator_builds = half_step_cache.base_operator_builds;
    result.full_step_ic0_factor_builds =
        full_step_cache.cached_cg.factor_builds();
    result.half_step_ic0_factor_builds =
        half_step_cache.cached_cg.factor_builds();
    result.temperature = std::move(state);
    result.final_time = time;
    return result;
}

} // namespace nova::fem
