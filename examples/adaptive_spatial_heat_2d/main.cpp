#include "nova/fem/error_indicators_2d.hpp"
#include "nova/fem/mesh_transfer_2d.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>

int main() {
    using namespace nova::fem;

    auto generated = make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 8U, 8U
    );

    TransientHeat2DProblem problem;
    problem.mesh = std::move(generated.mesh);
    problem.volumetric_capacity = [](double, double) {
        return 1.0;
    };
    problem.conductivity = [](double, double) {
        return 1.0;
    };
    problem.source = [](const double x, const double y, const double time) {
        const double dx = x - 0.5;
        const double dy = y - 0.5;
        return 40.0 * std::exp(-120.0 * (dx * dx + dy * dy)) *
            (1.0 + 0.1 * time);
    };
    problem.initial_temperature = [](double, double) {
        return 0.0;
    };
    for (const BoundaryEdge2D& edge : generated.boundary_edges) {
        problem.dirichlet_edges.push_back({
            edge,
            [](double, double, double) { return 0.0; }
        });
    }

    TransientHeat2DOptions options;
    options.order = H1TriangleOrder::P2;
    options.time_step = 0.05;
    options.steps = 1U;

    const auto coarse_solution = solve_transient_heat_2d(problem, options);
    const auto indicators = estimate_transient_heat_gradient_recovery_2d(
        problem, options.order, coarse_solution.temperature
    );
    const auto marked = mark_elements_dorfler_2d(
        indicators.element_errors, 0.5
    );
    if (marked.empty()) {
        throw std::runtime_error(
            "The localized thermal source produced no spatial refinement marks."
        );
    }

    // The accepted coarse solution becomes the exact P2 restart state.
    problem.initial_dof_values = coarse_solution.temperature;
    const auto refined = refine_transient_heat_problem_adaptively_2d(
        problem, options.order, marked
    );

    options.start_time = coarse_solution.final_time;
    const auto fine_solution = solve_transient_heat_2d(
        refined.problem, options
    );
    const double maximum_temperature = *std::max_element(
        fine_solution.temperature.begin(),
        fine_solution.temperature.end()
    );
    if (!std::isfinite(maximum_temperature)) {
        throw std::runtime_error("Adaptive thermal solution is non-finite.");
    }

    std::cout << "NOVA adaptive transient heat example\n"
              << "Coarse elements: " << problem.mesh.elements.size() << '\n'
              << "Marked elements: " << marked.size() << '\n'
              << "Refined elements: "
              << refined.problem.mesh.elements.size() << '\n'
              << "Recovery indicator: " << indicators.total_error << '\n'
              << "Maximum temperature: " << maximum_temperature << '\n'
              << "Final time: " << fine_solution.final_time << '\n';
    return 0;
}
