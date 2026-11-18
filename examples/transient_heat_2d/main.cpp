#include "nova/fem/transient_heat_2d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

int main() {
    // A warm square plate with its outer boundary held at zero.
    // Solve rho*c*dT/dt - div(k*grad(T)) = 0 for 0 <= t <= 0.1.
    auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 12U, 12U
    );
    nova::fem::TransientHeat2DProblem problem;
    problem.mesh = std::move(generated.mesh);
    problem.volumetric_capacity = [](double, double) { return 1.0; };
    problem.conductivity = [](double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    const double pi = std::acos(-1.0);
    problem.initial_temperature = [pi](double x, double y) {
        return std::sin(pi * x) * std::sin(pi * y);
    };
    for (const auto& edge : generated.boundary_edges) {
        problem.dirichlet_edges.push_back({
            edge, [](double, double, double) { return 0.0; }
        });
    }

    nova::fem::TransientHeat2DOptions options;
    options.order = nova::fem::H1TriangleOrder::P2;
    options.time_step = 0.01;
    options.steps = 10U;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;

    const auto result = nova::fem::solve_transient_heat_2d(problem, options);
    const auto extrema = std::minmax_element(
        result.temperature.begin(), result.temperature.end()
    );
    std::cout << std::fixed << std::setprecision(8)
              << "NOVA transient heat / Backward Euler\n"
              << "Order: P2\n"
              << "DOFs: " << result.temperature.size() << '\n'
              << "Time: " << result.final_time << '\n'
              << "Minimum temperature: " << *extrema.first << '\n'
              << "Maximum temperature: " << *extrema.second << '\n'
              << "Completed steps: " << result.linear_solvers.size() << '\n';
    return 0;
}
