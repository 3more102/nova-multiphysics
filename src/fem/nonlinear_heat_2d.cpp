#include "nova/fem/nonlinear_heat_2d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

// The linear thermal assembler requests k(x,y) at quadrature points.
// Locate the containing triangle and evaluate the current H1 FE field there.
// This intentionally simple O(number_of_elements) lookup per query is a
// reference implementation, not a production spatial search accelerator.
[[nodiscard]] double evaluate_fe_temperature(
    const H1FunctionSpace2D& space,
    const std::vector<double>& values,
    const double x,
    const double y
) {
    const TriangleMesh2D& mesh = space.mesh();
    constexpr double tolerance = 1.0e-11;
    for (std::size_t element_index = 0; element_index < mesh.elements.size();
         ++element_index) {
        const Triangle2D& cell = mesh.elements[element_index];
        const Point2D& a = mesh.nodes[cell[0]];
        const Point2D& b = mesh.nodes[cell[1]];
        const Point2D& c = mesh.nodes[cell[2]];
        const double j00 = b.x - a.x;
        const double j01 = c.x - a.x;
        const double j10 = b.y - a.y;
        const double j11 = c.y - a.y;
        const double determinant = j00 * j11 - j01 * j10;
        const double dx = x - a.x;
        const double dy = y - a.y;
        const double r = (j11 * dx - j01 * dy) / determinant;
        const double s = (j00 * dy - j10 * dx) / determinant;
        if (r < -tolerance || s < -tolerance ||
            r + s > 1.0 + tolerance) {
            continue;
        }
        const auto basis = evaluate_triangle_h1_basis(space.order(), r, s);
        const auto& dofs = space.dofs().element_dofs(element_index);
        double temperature = 0.0;
        for (std::size_t i = 0; i < basis.count; ++i) {
            temperature += basis.values[i] * values[dofs[i]];
        }
        if (!std::isfinite(temperature)) {
            throw std::overflow_error("Nonlinear H1 temperature evaluation overflowed.");
        }
        return temperature;
    }
    throw std::out_of_range(
        "Nonlinear thermal conductivity requested outside the FE mesh."
    );
}

void validate_options(
    const NonlinearHeat2DProblem& problem,
    const PicardHeat2DOptions& options
) {
    if (!problem.conductivity) {
        throw std::invalid_argument("Temperature-dependent conductivity is required.");
    }
    if (!std::isfinite(options.start_time) ||
        !std::isfinite(options.time_step) || options.time_step <= 0.0 ||
        options.steps == 0U) {
        throw std::invalid_argument("Invalid nonlinear thermal time window.");
    }
    if (!std::isfinite(options.absolute_tolerance) ||
        options.absolute_tolerance <= 0.0 ||
        !std::isfinite(options.relative_tolerance) ||
        options.relative_tolerance < 0.0) {
        throw std::invalid_argument("Invalid nonlinear Picard convergence tolerance.");
    }
    if (!std::isfinite(options.relaxation) ||
        options.relaxation <= 0.0 || options.relaxation > 1.0 ||
        options.max_iterations == 0U) {
        throw std::invalid_argument("Invalid nonlinear Picard iteration controls.");
    }
    const double end = options.start_time +
        options.time_step * static_cast<double>(options.steps);
    if (!std::isfinite(end) || end <= options.start_time) {
        throw std::invalid_argument("Nonlinear thermal end time is not representable.");
    }
}

[[nodiscard]] std::pair<double, double> update_error(
    const std::vector<double>& next,
    const std::vector<double>& guess,
    const PicardHeat2DOptions& options
) {
    if (next.size() != guess.size()) {
        throw std::runtime_error("Picard trial has an unexpected DOF count.");
    }
    double absolute = 0.0;
    double scaled = 0.0;
    for (std::size_t i = 0; i < next.size(); ++i) {
        if (!std::isfinite(next[i])) {
            throw std::runtime_error("Picard trial contains nonfinite DOFs.");
        }
        const double delta = std::abs(next[i] - guess[i]);
        const double scale = options.absolute_tolerance +
            options.relative_tolerance *
                std::max(std::abs(next[i]), std::abs(guess[i]));
        if (!std::isfinite(delta) || !std::isfinite(scale) || scale <= 0.0) {
            throw std::overflow_error("Picard convergence scaling overflowed.");
        }
        absolute = std::max(absolute, delta);
        scaled = std::max(scaled, delta / scale);
    }
    return {absolute, scaled};
}

} // namespace

PicardHeat2DResult solve_picard_heat_2d(
    const NonlinearHeat2DProblem& problem,
    const PicardHeat2DOptions& options
) {
    validate_options(problem, options);
    problem.thermal.mesh.validate();

    const H1FunctionSpace2D space(problem.thermal.mesh, options.order);
    std::vector<double> state = problem.thermal.initial_dof_values;
    if (state.empty()) {
        if (!problem.thermal.initial_temperature) {
            throw std::invalid_argument("Initial temperature or exact DOFs are required.");
        }
        state = space.interpolate(problem.thermal.initial_temperature);
    }
    if (state.size() != space.dimension() ||
        !std::all_of(state.begin(), state.end(),
            [](double value) { return std::isfinite(value); })) {
        throw std::invalid_argument("Invalid initial nonlinear thermal H1 DOFs.");
    }

    PicardHeat2DResult result;
    result.step_diagnostics.reserve(options.steps);

    TransientHeat2DOptions inner;
    inner.order = options.order;
    inner.scheme = TransientTimeScheme::BackwardEuler;
    inner.time_step = options.time_step;
    inner.steps = 1U;
    inner.solver_options = options.solver_options;

    for (std::size_t step = 0; step < options.steps; ++step) {
        const double current_time = options.start_time +
            options.time_step * static_cast<double>(step);
        inner.start_time = current_time;

        std::vector<double> guess = state;
        bool converged = false;
        PicardHeatStep2D diagnostic;
        diagnostic.time = current_time + options.time_step;

        for (std::size_t iteration = 0; iteration < options.max_iterations;
             ++iteration) {
            TransientHeat2DProblem frozen = problem.thermal;
            frozen.initial_dof_values = state;
            frozen.conductivity = [&](double x, double y) {
                const double t = evaluate_fe_temperature(space, guess, x, y);
                const double k = problem.conductivity(x, y, t);
                if (!std::isfinite(k) || k <= 0.0) {
                    throw std::invalid_argument(
                        "Nonlinear thermal conductivity must be finite and positive."
                    );
                }
                return k;
            };

            const auto trial = solve_transient_heat_2d(frozen, inner);
            if (trial.linear_solvers.size() != 1U) {
                throw std::runtime_error("Expected one linear thermal solve per Picard trial.");
            }
            result.linear_solvers.push_back(trial.linear_solvers.front());

            const auto [absolute, normalized] =
                update_error(trial.temperature, guess, options);
            diagnostic.iterations = iteration + 1U;
            diagnostic.max_absolute_update = absolute;
            diagnostic.normalized_update = normalized;

            if (normalized <= 1.0) {
                state = trial.temperature;
                converged = true;
                break;
            }

            for (std::size_t i = 0; i < guess.size(); ++i) {
                guess[i] += options.relaxation * (trial.temperature[i] - guess[i]);
                if (!std::isfinite(guess[i])) {
                    throw std::overflow_error("Relaxed Picard temperature is nonfinite.");
                }
            }
        }

        if (!converged) {
            throw std::runtime_error(
                "Nonlinear thermal Picard iterations failed at time step " +
                std::to_string(step + 1U) + "."
            );
        }
        result.step_diagnostics.push_back(diagnostic);
    }

    result.temperature = std::move(state);
    result.final_time = options.start_time +
        options.time_step * static_cast<double>(options.steps);
    return result;
}

} // namespace nova::fem
