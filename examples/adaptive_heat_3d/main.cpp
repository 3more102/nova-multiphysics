#include "nova/fem/adaptive_heat_3d.hpp"

#include "nova/fem/mesh3d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

int main() {
    auto grid = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U
    );

    nova::fem::TransientHeat3DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.volumetric_capacity = [](double, double, double) { return 1.0; };
    problem.conductivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double, double time) {
        return 1.0 + 0.25 * time;
    };
    problem.initial_temperature = [](double, double, double) { return 0.0; };
    for (const auto& face : grid.boundary_faces) {
        problem.dirichlet_faces.push_back({
            face, [](double, double, double, double) { return 0.0; }
        });
    }

    nova::fem::AdaptiveHeat3DOptions options;
    options.scheme = nova::fem::TransientTimeScheme::CrankNicolson;
    options.start_time = 0.0;
    options.end_time = 0.5;
    options.initial_step = 0.25;
    options.min_step = 1.0e-5;
    options.max_step = 0.25;
    options.absolute_tolerance = 1.0e-5;
    options.relative_tolerance = 1.0e-4;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;

    const auto solution = nova::fem::solve_adaptive_heat_3d(problem, options);
    const auto extrema = std::minmax_element(
        solution.temperature.begin(), solution.temperature.end()
    );
    double maximum_error = 0.0;
    double maximum_budget_residual = 0.0;
    for (const auto& attempt : solution.attempts) {
        if (attempt.accepted) {
            maximum_error = std::max(maximum_error, attempt.normalized_error);
        }
    }
    for (const auto& diagnostic : solution.diagnostics) {
        maximum_budget_residual = std::max(
            maximum_budget_residual, std::abs(diagnostic.balance_residual)
        );
    }

    std::cout << std::fixed << std::setprecision(8)
              << "NOVA P1 tetrahedral adaptive Crank-Nicolson heat example\n"
              << "Final time: " << solution.final_time << '\n'
              << "Tetrahedra: " << problem.mesh.elements.size() << '\n'
              << "P1 DOFs: " << solution.temperature.size() << '\n'
              << "Accepted steps: " << solution.accepted_steps << '\n'
              << "Rejected trials: " << solution.rejected_steps << '\n'
              << "Minimum temperature: " << *extrema.first << '\n'
              << "Maximum temperature: " << *extrema.second << '\n'
              << "Largest accepted normalized error: " << maximum_error << '\n'
              << "Maximum discrete energy-balance residual: "
              << maximum_budget_residual << '\n';
    return 0;
}
