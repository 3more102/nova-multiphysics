#include "nova/fem/adaptive_heat_3d.hpp"

#include "transient_heat_3d_internal.hpp"
#include "nova/linalg/conjugate_gradient.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {

void validate_options(const AdaptiveHeat3DOptions& options) {
    const double duration = options.end_time - options.start_time;
    if (!std::isfinite(options.start_time) || !std::isfinite(options.end_time) ||
        !std::isfinite(duration) || duration <= 0.0) {
        throw std::invalid_argument("Adaptive 3D heat needs finite increasing start/end times.");
    }
    if (!std::isfinite(options.initial_step) || options.initial_step <= 0.0 ||
        !std::isfinite(options.min_step) || options.min_step <= 0.0 ||
        !std::isfinite(options.max_step) || options.max_step < options.min_step) {
        throw std::invalid_argument("Adaptive 3D time-step bounds must be finite, positive and ordered.");
    }
    if (!std::isfinite(options.absolute_tolerance) ||
        options.absolute_tolerance <= 0.0 ||
        !std::isfinite(options.relative_tolerance) ||
        options.relative_tolerance < 0.0) {
        throw std::invalid_argument(
            "Adaptive 3D absolute tolerance must be positive and relative tolerance nonnegative."
        );
    }
    if (options.max_attempts == 0U) {
        throw std::invalid_argument("Adaptive 3D maximum attempts must be positive.");
    }
    if (options.scheme != TransientTimeScheme::BackwardEuler &&
        options.scheme != TransientTimeScheme::CrankNicolson) {
        throw std::invalid_argument("Adaptive 3D heat time scheme is unsupported.");
    }
}

[[nodiscard]] double normalized_step_error(
    const std::vector<double>& initial,
    const std::vector<double>& full,
    const std::vector<double>& half,
    const AdaptiveHeat3DOptions& options
) {
    if (initial.empty() || full.size() != initial.size() || half.size() != initial.size()) {
        throw std::invalid_argument("Adaptive 3D trial field dimension mismatch.");
    }
    const double divisor =
        options.scheme == TransientTimeScheme::BackwardEuler ? 1.0 : 3.0;
    double sum = 0.0;
    for (std::size_t i = 0; i < initial.size(); ++i) {
        const double scale = options.absolute_tolerance +
            options.relative_tolerance * std::max(std::abs(initial[i]), std::abs(half[i]));
        if (!std::isfinite(scale) || scale <= 0.0) {
            throw std::overflow_error("Adaptive 3D error scale is non-finite or non-positive.");
        }
        const double error = ((half[i] - full[i]) / divisor) / scale;
        if (!std::isfinite(error)) {
            throw std::overflow_error("Adaptive 3D local truncation error is non-finite.");
        }
        sum += error * error;
        if (!std::isfinite(sum)) {
            throw std::overflow_error("Adaptive 3D squared error accumulation overflowed.");
        }
    }
    return std::sqrt(sum / static_cast<double>(initial.size()));
}

[[nodiscard]] TransientHeatStepDiagnostics3D combine_half_budgets(
    const std::vector<TransientHeatStepDiagnostics3D>& half
) {
    if (half.size() != 2U) {
        throw std::logic_error("Adaptive 3D heat expected exactly two half-step budgets.");
    }
    auto diagnostic = half.back();
    diagnostic.source_energy += half.front().source_energy;
    diagnostic.neumann_energy += half.front().neumann_energy;
    diagnostic.robin_energy += half.front().robin_energy;
    diagnostic.dirichlet_energy += half.front().dirichlet_energy;
    diagnostic.balance_residual += half.front().balance_residual;
    return diagnostic;
}

} // namespace

AdaptiveHeat3DResult solve_adaptive_heat_3d(
    const TransientHeat3DProblem& problem,
    const AdaptiveHeat3DOptions& options
) {
    validate_options(options);

    // Topology, boundary assignments, material operators and element Jacobians
    // are prepared once for the entire adaptive study, including rejected trials.
    const auto prepared = detail::prepare_transient_heat_3d(problem);
    std::vector<double> state = detail::initial_transient_heat_state_3d(problem);

    // Separate slots avoid repeatedly evicting the full-step factor when
    // computing two half steps. Exact matrix equality controls invalidation.
    linalg::CachedConjugateGradient full_step_cache;
    linalg::CachedConjugateGradient half_step_cache;
    AdaptiveHeat3DResult result;
    double time = options.start_time;
    double step = std::clamp(options.initial_step, options.min_step, options.max_step);
    const double method_order =
        options.scheme == TransientTimeScheme::BackwardEuler ? 1.0 : 2.0;

    while (time < options.end_time) {
        if (result.attempts.size() >= options.max_attempts) {
            throw std::runtime_error("Adaptive 3D heat maximum step attempts exhausted.");
        }
        const double remaining = options.end_time - time;
        const double dt = std::min(step, remaining);
        if (!std::isfinite(dt) || dt <= 0.0 || time + dt <= time) {
            throw std::runtime_error(
                "Adaptive 3D time step is too small to advance floating-point time."
            );
        }

        TransientHeat3DOptions fixed;
        fixed.scheme = options.scheme;
        fixed.start_time = time;
        fixed.time_step = dt;
        fixed.steps = 1U;
        fixed.solver_options = options.solver_options;
        const auto full = detail::solve_prepared_transient_heat_3d(
            problem, fixed, prepared, state, &full_step_cache
        );

        fixed.time_step = 0.5 * dt;
        fixed.steps = 2U;
        auto half = detail::solve_prepared_transient_heat_3d(
            problem, fixed, prepared, state, &half_step_cache
        );

        const double error = normalized_step_error(
            state, full.temperature, half.temperature, options
        );
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
        } else {
            ++result.rejected_steps;
            if (dt <= options.min_step) {
                throw std::runtime_error(
                    "Adaptive 3D heat error exceeds tolerance at minimum time step."
                );
            }
        }

        const double bounded_error = std::max(error, 1.0e-14);
        double factor = 0.9 * std::pow(bounded_error, -1.0 / (method_order + 1.0));
        factor = std::clamp(factor, 0.2, accepted ? 4.0 : 0.9);
        step = std::clamp(dt * factor, options.min_step, options.max_step);
    }

    result.temperature = std::move(state);
    result.full_step_ic0_factor_builds = full_step_cache.factor_builds();
    result.half_step_ic0_factor_builds = half_step_cache.factor_builds();
    result.final_time = time;
    return result;
}

} // namespace nova::fem
