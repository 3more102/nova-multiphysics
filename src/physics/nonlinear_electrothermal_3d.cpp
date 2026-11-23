#include "nova/physics/nonlinear_electrothermal_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace nova::physics {
namespace {

void validate_options(const NonlinearElectrothermal3DProblem& problem,
                      const NonlinearElectrothermal3DOptions& options) {
    if (!problem.conductivity_scale) {
        throw std::invalid_argument("Electrothermal feedback law is required.");
    }
    if (!std::isfinite(problem.initial_temperature_kelvin) ||
        problem.initial_temperature_kelvin <= 0.0) {
        throw std::invalid_argument("Electrothermal initial temperature must be positive and finite.");
    }
    if (options.max_iterations == 0 ||
        !std::isfinite(options.relaxation) ||
        options.relaxation <= 0.0 || options.relaxation > 1.0 ||
        !std::isfinite(options.absolute_tolerance_kelvin) ||
        options.absolute_tolerance_kelvin < 0.0 ||
        !std::isfinite(options.relative_tolerance) ||
        options.relative_tolerance < 0.0 ||
        (options.absolute_tolerance_kelvin == 0.0 &&
         options.relative_tolerance == 0.0) ||
        (options.adaptive_relaxation &&
         (!std::isfinite(options.minimum_relaxation) ||
          options.minimum_relaxation <= 0.0 ||
          options.minimum_relaxation > options.relaxation))) {
        throw std::invalid_argument("Invalid electrothermal Picard solver options.");
    }

    const auto require_cellwise_material = [](const fem::ScalarDiffusion3DProblem& p) {
        const auto cells = p.mesh.elements.size();
        const bool scalar = !p.cell_diffusivity.empty();
        const bool tensor = !p.cell_diffusivity_tensors.empty();
        if (scalar == tensor || p.diffusivity || p.diffusivity_tensor ||
            (scalar && p.cell_diffusivity.size() != cells) ||
            (tensor && p.cell_diffusivity_tensors.size() != cells)) {
            throw std::invalid_argument(
                "Temperature-dependent conductivity requires exactly one "
                "complete per-cell scalar or tensor reference material.");
        }
    };
    require_cellwise_material(problem.base.electrical);
    // Preserve legacy spatial/tensor thermal material paths unless k(T)
    // is explicitly requested.
    if (problem.thermal_conductivity_scale) {
        require_cellwise_material(problem.base.thermal);
    }
}

double checked_product(double base, double scale) {
    const double value = base * scale;
    if (!std::isfinite(value)) {
        throw std::overflow_error("Electrothermal feedback conductivity overflows.");
    }
    return value;
}

void update_scaled_cell_material(
    const fem::ScalarDiffusion3DProblem& original,
    fem::ScalarDiffusion3DProblem& updated,
    const std::vector<double>& temperature,
    const std::function<double(double)>& scaling_law
) {
    for (std::size_t i = 0; i < original.mesh.elements.size(); ++i) {
        const auto& tetrahedron = original.mesh.elements[i];
        long double sum = 0.0L;
        for (const auto vertex : tetrahedron) {
            if (vertex >= temperature.size()) {
                throw std::invalid_argument(
                    "Electrothermal material mesh references an invalid vertex.");
            }
            sum += static_cast<long double>(temperature[vertex]);
        }
        const double mean_temperature = static_cast<double>(sum / 4.0L);
        if (!std::isfinite(mean_temperature)) {
            throw std::overflow_error("Electrothermal cell temperature is not finite.");
        }
        const double scale = scaling_law(mean_temperature);
        if (!std::isfinite(scale) || scale <= 0.0) {
            throw std::invalid_argument(
                "Electrothermal conductivity multiplier must be positive and finite.");
        }
        if (!original.cell_diffusivity.empty()) {
            const double coefficient = checked_product(original.cell_diffusivity[i], scale);
            if (coefficient <= 0.0) {
                throw std::invalid_argument(
                    "Electrothermal scaled scalar conductivity must be positive.");
            }
            updated.cell_diffusivity[i] = coefficient;
        } else {
            const auto& a = original.cell_diffusivity_tensors[i];
            auto& b = updated.cell_diffusivity_tensors[i];
            b = {
                checked_product(a.xx, scale), checked_product(a.yy, scale),
                checked_product(a.zz, scale), checked_product(a.xy, scale),
                checked_product(a.xz, scale), checked_product(a.yz, scale)
            };
            // Shared stationary FEM assembly validates tensor SPD.
        }
    }
}

} // namespace

NonlinearElectrothermal3DResult solve_nonlinear_electrothermal_3d(
    const NonlinearElectrothermal3DProblem& problem,
    const NonlinearElectrothermal3DOptions& options,
    const linalg::SolverOptions& electrical_options,
    const linalg::SolverOptions& thermal_options
) {
    validate_options(problem, options);
    const auto n = problem.base.electrical.mesh.nodes.size();
    std::vector<double> temperature(n, problem.initial_temperature_kelvin);

    NonlinearElectrothermal3DResult result;
    result.history.reserve(options.max_iterations);
    double relaxation = options.relaxation;
    for (std::size_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
        auto trial = problem.base;
        update_scaled_cell_material(
            problem.base.electrical, trial.electrical, temperature,
            problem.conductivity_scale);
        if (problem.thermal_conductivity_scale) {
            update_scaled_cell_material(
                problem.base.thermal, trial.thermal, temperature,
                problem.thermal_conductivity_scale);
        }
        result.fields = solve_one_way_electrothermal_3d(
            trial, electrical_options, thermal_options);
        const auto& solved = result.fields.thermal_solution.values;
        if (solved.size() != temperature.size()) {
            throw std::runtime_error("Electrothermal temperature field size mismatch.");
        }
        result.iterations = iteration;
        result.last_temperature_iterate = temperature;
        result.max_fixed_point_residual_kelvin = 0.0;
        result.max_normalized_residual = 0.0;
        for (std::size_t node = 0; node < n; ++node) {
            if (!std::isfinite(solved[node])) {
                throw std::runtime_error("Electrothermal solved temperature is not finite.");
            }
            const double diff = std::abs(solved[node] - temperature[node]);
            const double scale = std::max(std::abs(solved[node]), std::abs(temperature[node]));
            const double threshold = options.absolute_tolerance_kelvin +
                                     options.relative_tolerance * scale;
            if (!std::isfinite(diff) || !std::isfinite(threshold)) {
                throw std::overflow_error("Electrothermal convergence metric overflows.");
            }
            result.max_fixed_point_residual_kelvin =
                std::max(result.max_fixed_point_residual_kelvin, diff);
            result.max_normalized_residual = std::max(
                result.max_normalized_residual, diff / threshold);
        }
        if (!std::isfinite(result.max_normalized_residual)) {
            throw std::overflow_error("Electrothermal normalized residual overflows.");
        }
        if (solved.empty()) {
            throw std::runtime_error("Electrothermal temperature field is empty.");
        }
        const auto [minimum, maximum] = std::minmax_element(
            solved.begin(), solved.end());
        NonlinearElectrothermal3DIteration diagnostic;
        diagnostic.iteration = iteration;
        diagnostic.max_fixed_point_residual_kelvin =
            result.max_fixed_point_residual_kelvin;
        diagnostic.max_normalized_residual = result.max_normalized_residual;
        diagnostic.total_joule_power_watts = result.fields.total_joule_power_watts;
        diagnostic.signed_thermal_balance_defect_watts =
            result.fields.thermal_balance.signed_balance_defect;
        diagnostic.min_solved_temperature_kelvin = *minimum;
        diagnostic.max_solved_temperature_kelvin = *maximum;
        if (result.max_normalized_residual <= 1.0) {
            result.converged = true;
            result.history.push_back(diagnostic);
            return result;
        }
        if (iteration == options.max_iterations) {
            result.history.push_back(diagnostic);
            break;
        }
        if (options.adaptive_relaxation && !result.history.empty()) {
            const long double current =
                static_cast<long double>(result.max_normalized_residual);
            const long double previous =
                static_cast<long double>(result.history.back().max_normalized_residual);
            if (current > 1.05L * previous) {
                relaxation = std::max(options.minimum_relaxation, 0.5 * relaxation);
            } else if (current < 0.8L * previous) {
                relaxation = std::min(options.relaxation, 1.2 * relaxation);
            }
        }
        diagnostic.relaxation_applied = relaxation;
        result.history.push_back(diagnostic);
        for (std::size_t node = 0; node < n; ++node) {
            // Convex update limits overshoot; no state from the caller is changed.
            temperature[node] += relaxation * (solved[node] - temperature[node]);
            if (!std::isfinite(temperature[node])) {
                throw std::overflow_error("Electrothermal relaxed temperature overflows.");
            }
        }
    }
    return result;
}

} // namespace nova::physics
