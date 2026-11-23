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
         options.relative_tolerance == 0.0)) {
        throw std::invalid_argument("Invalid electrothermal Picard solver options.");
    }

    const auto& e = problem.base.electrical;
    const auto cells = e.mesh.elements.size();
    const bool scalar = !e.cell_diffusivity.empty();
    const bool tensor = !e.cell_diffusivity_tensors.empty();
    if (scalar == tensor || e.diffusivity || e.diffusivity_tensor ||
        (scalar && e.cell_diffusivity.size() != cells) ||
        (tensor && e.cell_diffusivity_tensors.size() != cells)) {
        throw std::invalid_argument(
            "Electrothermal feedback requires exactly one complete cellwise "
            "scalar or tensor reference conductivity field.");
    }
}

double checked_product(double base, double scale) {
    const double value = base * scale;
    if (!std::isfinite(value)) {
        throw std::overflow_error("Electrothermal feedback conductivity overflows.");
    }
    return value;
}

void update_electrical_material(
    const NonlinearElectrothermal3DProblem& problem,
    const std::vector<double>& temperature,
    fem::ScalarDiffusion3DProblem& electrical
) {
    const auto& original = problem.base.electrical;
    for (std::size_t i = 0; i < original.mesh.elements.size(); ++i) {
        const auto& t = original.mesh.elements[i];
        long double sum = 0.0L;
        for (const auto vertex : t) {
            if (vertex >= temperature.size()) {
                throw std::invalid_argument("Electrothermal mesh references an invalid vertex.");
            }
            sum += static_cast<long double>(temperature[vertex]);
        }
        const double cell_temperature = static_cast<double>(sum / 4.0L);
        if (!std::isfinite(cell_temperature)) {
            throw std::overflow_error("Electrothermal cell temperature is not finite.");
        }
        const double scale = problem.conductivity_scale(cell_temperature);
        if (!std::isfinite(scale) || scale <= 0.0) {
            throw std::invalid_argument(
                "Electrothermal conductivity scale must be positive and finite.");
        }
        if (!original.cell_diffusivity.empty()) {
            electrical.cell_diffusivity[i] =
                checked_product(original.cell_diffusivity[i], scale);
            if (electrical.cell_diffusivity[i] <= 0.0) {
                throw std::invalid_argument(
                    "Electrothermal scaled scalar conductivity must be positive.");
            }
        } else {
            const auto& a = original.cell_diffusivity_tensors[i];
            auto& b = electrical.cell_diffusivity_tensors[i];
            b = {
                checked_product(a.xx, scale), checked_product(a.yy, scale),
                checked_product(a.zz, scale), checked_product(a.xy, scale),
                checked_product(a.xz, scale), checked_product(a.yz, scale)
            };
            // SPD validation remains the responsibility of the shared
            // stationary FEM conductivity assembly at quadrature points.
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
    for (std::size_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
        auto trial = problem.base;
        update_electrical_material(problem, temperature, trial.electrical);
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
        if (result.max_normalized_residual <= 1.0) {
            result.converged = true;
            return result;
        }
        if (iteration == options.max_iterations) break;
        for (std::size_t node = 0; node < n; ++node) {
            // Convex update limits overshoot; no state from the caller is changed.
            temperature[node] += options.relaxation * (solved[node] - temperature[node]);
            if (!std::isfinite(temperature[node])) {
                throw std::overflow_error("Electrothermal relaxed temperature overflows.");
            }
        }
    }
    return result;
}

} // namespace nova::physics
