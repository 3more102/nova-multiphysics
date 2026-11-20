#include "nova/fem/adaptive_heat_3d.hpp"

#include "transient_heat_3d_internal.hpp"

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

    // Static mesh/material/boundary integration is performed once for the
    // complete adaptive study, including all rejected trials.
    const auto prepared = detail::prepare_transient_heat_3d(problem);
    std::vector<double> state = prepared.initial_state;

    // Full and half steps have different dt-dependent matrices. Separate
    // slots retain each operator set and its exact-matrix IC(0) factor.
    detail::TransientHeat3DStepCache full_step_cache;
    detail::TransientHeat3DStepCache half_step_cache;

    AdaptiveHeat3DResult result;
    double time = options.start_time;
    double step = std::clamp(
        options.initial_step, options.min_step, options.max_step);
    const double order =
        options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 2.0;

    if (options.on_step) {
        options.on_step(0U, time, state);
    }

    while (time < options.end_time) {
        if (result.attempts.size() >= options.max_attempts) {
            throw std::runtime_error("Adaptive 3D maximum attempts exceeded.");
        }
        const double dt = std::min(step, options.end_time - time);
        if (!(dt > 0.0) || !std::isfinite(dt) ||
            !std::isfinite(time + dt) || !(time + dt > time)) {
            throw std::runtime_error(
                "Adaptive 3D time cannot advance representably.");
        }

        TransientHeat3DOptions fixed;
        fixed.scheme = options.scheme;
        fixed.start_time = time;
        fixed.time_step = dt;
        fixed.steps = 1U;
        fixed.solver_options = options.solver_options;
        fixed.retain_step_solutions = false;
        const auto full = detail::solve_prepared_transient_heat_3d(
            prepared, fixed, state, full_step_cache);

        fixed.time_step = 0.5 * dt;
        fixed.steps = 2U;
        auto half = detail::solve_prepared_transient_heat_3d(
            prepared, fixed, state, half_step_cache);

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

    result.full_step_operator_builds = full_step_cache.operator_builds;
    result.half_step_operator_builds = half_step_cache.operator_builds;
    result.full_step_ic0_factor_builds =
        full_step_cache.cached_cg.factor_builds();
    result.half_step_ic0_factor_builds =
        half_step_cache.cached_cg.factor_builds();
    result.temperature = std::move(state);
    result.final_time = time;
    return result;
}

} // namespace nova::fem
