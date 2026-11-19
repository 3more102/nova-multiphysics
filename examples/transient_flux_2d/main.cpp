#include "nova/fem/transient_heat_2d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

int main() {
    // Heat enters the left side as a prescribed inward conductive flux.
    // The right side exchanges heat with a warmer environment by convection.
    // The horizontal sides are insulated.
    auto grid = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 10U, 10U
    );
    nova::fem::TransientHeat2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.volumetric_capacity = [](double, double) { return 2.0; };
    problem.conductivity = [](double, double) { return 0.5; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.initial_temperature = [](double, double) { return 0.0; };

    for (const auto& edge : grid.boundary_edges) {
        const auto& a = problem.mesh.nodes[edge[0]];
        const auto& b = problem.mesh.nodes[edge[1]];
        if (a.x == 0.0 && b.x == 0.0) {
            problem.neumann_edges.push_back({
                edge,
                [](double, double, double) { return -0.2; } // inward heating
            });
        }
        if (a.x == 1.0 && b.x == 1.0) {
            problem.robin_edges.push_back({
                edge,
                [](double, double, double) { return 4.0; },
                [](double, double, double t) { return 1.0 + 0.5 * t; }
            });
        }
    }

    nova::fem::TransientHeat2DOptions options;
    options.order = nova::fem::H1TriangleOrder::P2;
    options.scheme = nova::fem::TransientTimeScheme::CrankNicolson;
    options.steps = 10U;
    options.time_step = 0.02;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;

    const auto result = nova::fem::solve_transient_heat_2d(problem, options);
    const auto extrema = std::minmax_element(
        result.temperature.begin(), result.temperature.end()
    );
    double boundary_energy = 0.0;
    double max_balance_error = 0.0;
    for (const auto& step : result.diagnostics) {
        boundary_energy += step.neumann_energy + step.robin_energy;
        max_balance_error = std::max(
            max_balance_error, std::abs(step.balance_residual)
        );
    }

    std::cout << std::fixed << std::setprecision(9)
              << "NOVA transient boundary-flux example\n"
              << "Method: P2 Crank-Nicolson\n"
              << "Final time: " << result.final_time << '\n'
              << "Degrees of freedom: " << result.temperature.size() << '\n'
              << "Min temperature: " << *extrema.first << '\n'
              << "Max temperature: " << *extrema.second << '\n'
              << "Final heat content: " << result.diagnostics.back().heat_content << '\n'
              << "Integrated boundary heat input: " << boundary_energy << '\n'
              << "Max discrete energy-balance residual: " << max_balance_error << '\n';
    return 0;
}
