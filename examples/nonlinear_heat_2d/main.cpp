#include "nova/fem/nonlinear_heat_2d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <utility>

int main() {
    // Manufactured nonlinear transient field:
    // T(x,y,t)=x+y+t, k(T)=1+0.5*T^2, rho*c=1.
    // Q = rho*c*T_t - div(k grad(T)) = 1 - 2*T.
    auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 6U, 6U
    );
    nova::fem::NonlinearHeat2DProblem problem;
    problem.mesh = std::move(generated.mesh);
    problem.volumetric_capacity = [](double, double) { return 1.0; };
    problem.conductivity = [](double, double, double t) {
        return 1.0 + 0.5 * t * t;
    };
    problem.conductivity_temperature_derivative =
        [](double, double, double t) { return t; };
    problem.source = [](double x, double y, double t) {
        return 1.0 - 2.0 * (x + y + t);
    };
    problem.initial_temperature = [](double x, double y) { return x + y; };

    for (const auto& edge : generated.boundary_edges) {
        problem.dirichlet_edges.push_back({
            edge, [](double x, double y, double t) {
                return x + y + t;
            }
        });
    }

    nova::fem::NonlinearHeat2DOptions options;
    options.order = nova::fem::H1TriangleOrder::P2;
    options.steps = 5U;
    options.time_step = 0.1;
    options.newton_options.gmres.solver.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;

    const auto result = nova::fem::solve_nonlinear_heat_2d(
        problem, options
    );
    nova::fem::H1FunctionSpace2D space(problem.mesh, options.order);
    const auto exact = space.interpolate([&](double x, double y) {
        return x + y + result.final_time;
    });
    double max_error = 0.0;
    std::size_t total_newton_iterations = 0U;
    for (std::size_t i = 0; i < exact.size(); ++i) {
        max_error = std::max(
            max_error, std::abs(result.temperature[i] - exact[i])
        );
    }
    for (const auto& step : result.newton_steps) {
        total_newton_iterations += step.iterations;
    }

    std::cout << std::fixed << std::setprecision(10)
              << "NOVA nonlinear transient thermal P2 verification\n"
              << "Final time: " << result.final_time << '\n'
              << "H1 degrees of freedom: " << result.temperature.size() << '\n'
              << "Newton time steps: " << result.newton_steps.size() << '\n'
              << "Total Newton iterations: " << total_newton_iterations << '\n'
              << "Maximum analytical DOF error: " << max_error << '\n';
    return max_error < 1.0e-7 ? 0 : 1;
}
