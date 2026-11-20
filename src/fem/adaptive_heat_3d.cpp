#include "nova/fem/adaptive_heat_3d.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "transient_heat_3d_internal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {

void validate(const AdaptiveHeat3DOptions& options) {
    const double duration = options.end_time - options.start_time;
    if (!std::isfinite(options.start_time) || !std::isfinite(options.end_time) ||
        !std::isfinite(duration) || duration <= 0.0) {
        throw std::invalid_argument("Adaptive 3D heat needs finite increasing times.");
    }
    if (!std::isfinite(options.initial_step) || options.initial_step <= 0.0 ||
        !std::isfinite(options.min_step) || options.min_step <= 0.0 ||
        !std::isfinite(options.max_step) || options.max_step < options.min_step) {
        throw std::invalid_argument("Adaptive 3D heat step bounds must be finite, positive and ordered.");
    }
    if (!std::isfinite(options.absolute_tolerance) || options.absolute_tolerance < 0.0 ||
        !std::isfinite(options.relative_tolerance) || options.relative_tolerance < 0.0 ||
        (options.absolute_tolerance == 0.0 && options.relative_tolerance == 0.0)) {
        throw std::invalid_argument("Adaptive 3D heat needs a positive error tolerance.");
    }
    if (options.max_attempts == 0U) {
        throw std::invalid_argument("Adaptive 3D heat requires at least one attempt.");
    }
    if (options.scheme != TransientHeat3DScheme::BackwardEuler &&
        options.scheme != TransientHeat3DScheme::CrankNicolson) {
        throw std::invalid_argument("Unsupported adaptive 3D time scheme.");
    }
}

double error_ratio(const std::vector<double>& old,
                   const std::vector<double>& full,
                   const std::vector<double>& fine,
                   const AdaptiveHeat3DOptions& options) {
    if (old.empty() || full.size() != old.size() || fine.size() != old.size()) {
        throw std::logic_error("Adaptive 3D heat trial state size mismatch.");
    }
    const double divisor =
        options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 3.0;
    double worst = 0.0;
    for (std::size_t i = 0; i < old.size(); ++i) {
        const double scale = options.absolute_tolerance +
            options.relative_tolerance *
                std::max({std::abs(old[i]), std::abs(full[i]), std::abs(fine[i])});
        if (!(scale > 0.0) || !std::isfinite(scale)) {
            throw std::overflow_error("Adaptive 3D heat error scale is invalid.");
        }
        const double normalized = ((fine[i] - full[i]) / divisor) / scale;
        if (!std::isfinite(normalized)) {
            throw std::overflow_error("Adaptive 3D heat error estimate overflowed.");
        }
        worst = std::max(worst, std::abs(normalized));
    }
    return worst;
}

TransientHeatStepDiagnostics3D combine_half_budgets(
    const TransientHeatStepDiagnostics3D& first,
    const TransientHeatStepDiagnostics3D& second
) {
    auto combined = second;
    combined.source_energy += first.source_energy;
    combined.neumann_energy += first.neumann_energy;
    combined.robin_energy += first.robin_energy;
    combined.dirichlet_energy += first.dirichlet_energy;
    combined.balance_residual += first.balance_residual;
    if (!std::isfinite(combined.source_energy) ||
        !std::isfinite(combined.neumann_energy) ||
        !std::isfinite(combined.robin_energy) ||
        !std::isfinite(combined.dirichlet_energy) ||
        !std::isfinite(combined.balance_residual)) {
        throw std::overflow_error("Adaptive 3D heat combined budget overflowed.");
    }
    return combined;
}

} // namespace

AdaptiveHeat3DResult solve_adaptive_heat_3d(
    const TransientHeat3DProblem& problem, const AdaptiveHeat3DOptions& options
) {
    validate(options);
    AdaptiveHeat3DResult result;
    nova::linalg::CachedConjugateGradient full_step_cache;
    nova::linalg::CachedConjugateGradient half_step_cache;
    std::vector<double> state = problem.initial_temperature;
    double time = options.start_time;
    double proposed = std::clamp(options.initial_step, options.min_step, options.max_step);
    const double order = options.scheme == TransientHeat3DScheme::BackwardEuler ? 1.0 : 2.0;
    bool initial_emitted = false;

    while (time < options.end_time) {
        if (result.attempts.size() >= options.max_attempts) {
            throw std::runtime_error("Adaptive 3D heat maximum attempts exhausted.");
        }
        const double remaining = options.end_time - time;
        const double dt = std::min(proposed, remaining);
        const double half_dt = 0.5 * dt;
        const double midpoint = time + half_dt;
        const double target = time + dt;
        if (!std::isfinite(dt) || !(half_dt > 0.0) ||
            !(midpoint > time) || !(target > midpoint) || !std::isfinite(target)) {
            throw std::runtime_error("Adaptive 3D heat cannot advance representable time.");
        }

        TransientHeat3DProblem trial = problem;
        trial.initial_temperature = state;
        TransientHeat3DOptions fixed;
        fixed.scheme = options.scheme;
        fixed.start_time = time;
        fixed.time_step = dt;
        fixed.steps = 1U;
        fixed.solver_options = options.solver_options;
        fixed.retain_step_solutions = false;

        if (!initial_emitted) {
            fixed.on_step = [&](std::size_t index, double physical_time,
                                const std::vector<double>& projected) {
                if (index == 0U) {
                    state = projected; // exact Dirichlet-projected restart DOFs
                    initial_emitted = true;
                    if (options.on_step) options.on_step(0U, physical_time, state);
                }
            };
        }
        const auto full = solve_transient_heat_3d_with_cache(trial, fixed, full_step_cache);
        fixed.on_step = {};

        // Evaluate each half step from its own exact accepted candidate,
        // not from a CG defect-correction vector.
        trial.initial_temperature = state;
        fixed.time_step = half_dt;
        auto first_half = solve_transient_heat_3d_with_cache(trial, fixed, half_step_cache);
        trial.initial_temperature = first_half.temperature;
        fixed.start_time = first_half.final_time;
        auto second_half = solve_transient_heat_3d_with_cache(trial, fixed, half_step_cache);

        const double error = error_ratio(
            state, full.temperature, second_half.temperature, options);
        const bool accepted = error <= 1.0;
        result.attempts.push_back({time, dt, error, accepted});
        if (accepted) {
            state = second_half.temperature;
            result.diagnostics.push_back(combine_half_budgets(
                first_half.diagnostics.at(0), second_half.diagnostics.at(0)));
            for (auto& solver : first_half.linear_solvers) {
                result.substep_solvers.push_back(std::move(solver));
            }
            for (auto& solver : second_half.linear_solvers) {
                result.substep_solvers.push_back(std::move(solver));
            }
            time = dt == remaining ? options.end_time : target;
            ++result.accepted_steps;
            if (options.on_step) options.on_step(result.accepted_steps, time, state);
        } else {
            ++result.rejected_steps;
            if (dt <= options.min_step) {
                throw std::runtime_error(
                    "Adaptive 3D heat tolerance cannot be met at minimum step.");
            }
        }

        const double safe_error = std::max(error, 1.0e-14);
        double factor = 0.9 * std::pow(safe_error, -1.0 / (order + 1.0));
        factor = std::clamp(factor, 0.2, accepted ? 2.0 : 0.9);
        proposed = std::clamp(dt * factor, options.min_step, options.max_step);
    }
    result.full_step_ic0_factor_builds = full_step_cache.factor_builds();
    result.half_step_ic0_factor_builds = half_step_cache.factor_builds();
    result.final_time = time;
    result.temperature = std::move(state);
    return result;
}

} // namespace nova::fem
