#include "nova/fem/adaptive_heat_2d.hpp"
#include "nova/io/vtk_legacy_2d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <utility>

int main(const int argc, char* argv[]) {
    if (argc != 1 && (argc != 3 ||
        std::string_view(argv[1]) != "--vtk")) {
        std::cerr << "Usage: nova_adaptive_heat_2d [--vtk output.vtk]\\n";
        return 2;
    }
    auto grid = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 6U, 6U
    );

    nova::fem::TransientHeat2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.volumetric_capacity = [](double, double) { return 1.0; };
    problem.conductivity = [](double, double) { return 1.0; };
    problem.source = [](double, double, double t) {
        return 1.0 + 3.0 * t * t;
    };
    problem.initial_temperature = [](double, double) { return 0.0; };
    for (const auto& edge : grid.boundary_edges) {
        problem.dirichlet_edges.push_back({
            edge,
            [](double, double, double) { return 0.0; }
        });
    }

    nova::fem::AdaptiveHeat2DOptions options;
    options.order = nova::fem::H1TriangleOrder::P2;
    options.scheme = nova::fem::TransientTimeScheme::CrankNicolson;
    options.start_time = 0.0;
    options.end_time = 1.0;
    options.initial_step = 0.4;
    options.min_step = 1.0e-5;
    options.max_step = 0.4;
    options.absolute_tolerance = 1.0e-5;
    options.relative_tolerance = 1.0e-4;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;

    const auto solution = nova::fem::solve_adaptive_heat_2d(
        problem, options
    );
    const auto extrema = std::minmax_element(
        solution.temperature.begin(), solution.temperature.end()
    );
    double max_normalized_error = 0.0;
    double max_budget_residual = 0.0;
    for (const auto& attempt : solution.attempts) {
        if (attempt.accepted) {
            max_normalized_error = std::max(
                max_normalized_error, attempt.normalized_error
            );
        }
    }
    for (const auto& diagnostic : solution.diagnostics) {
        max_budget_residual = std::max(
            max_budget_residual, std::abs(diagnostic.balance_residual)
        );
    }

    std::cout << std::fixed << std::setprecision(8)
              << "NOVA P2 adaptive Crank-Nicolson heat example\n"
              << "Final time: " << solution.final_time << '\n'
              << "H1 DOFs: " << solution.temperature.size() << '\n'
              << "Accepted steps: " << solution.accepted_steps << '\n'
              << "Rejected trials: " << solution.rejected_steps << '\n'
              << "Minimum temperature: " << *extrema.first << '\n'
              << "Maximum temperature: " << *extrema.second << '\n'
              << "Largest accepted normalized error: " << max_normalized_error << '\n'
              << "Maximum discrete energy-balance residual: "
              << max_budget_residual << '\n';
    if (argc == 3) {
        std::ofstream vtk_output(argv[2]);
        if (!vtk_output) {
            std::cerr << "Cannot open VTK output file: " << argv[2] << '\\n';
            return 2;
        }
        const nova::fem::H1FunctionSpace2D space(
            problem.mesh, options.order
        );
        nova::io::write_legacy_vtk_scalar_2d(
            vtk_output, space, solution.temperature
        );
        std::cout << "VTK result saved: " << argv[2] << '\\n';
    }
    return 0;
}
