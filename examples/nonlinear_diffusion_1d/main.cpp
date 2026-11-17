#include "nova/fem/nonlinear_diffusion_1d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {

constexpr double pi = 3.14159265358979323846;

double exact_solution(const double x) {
    return std::sin(pi * x);
}

double manufactured_source(const double x) {
    const double u = exact_solution(x);
    const double derivative = pi * std::cos(pi * x);
    const double second_derivative = -pi * pi * u;

    // f = -d/dx[(1+u^2)u'].
    return -(
        2.0 * u * derivative * derivative +
        (1.0 + u * u) * second_derivative
    );
}

} // namespace

int main() {
    nova::fem::NonlinearDiffusion1DProblem problem;
    problem.length = 1.0;
    problem.elements = 32;
    problem.source = manufactured_source;
    problem.left_value = 0.0;
    problem.right_value = 0.0;

    nova::nonlinear::NewtonOptions options;
    options.residual_tolerance = 1.0e-10;
    options.gmres.solver.relative_tolerance = 1.0e-11;
    options.gmres.solver.absolute_tolerance = 1.0e-13;

    const auto result = nova::fem::solve_nonlinear_diffusion_1d(
        problem,
        {},
        options
    );

    std::cout << "NOVA - Nonlinear 1D diffusion verification\n";
    std::cout << "Equation: -d/dx[(1+u^2)u'] = f(x), x in [0,1]\n";
    std::cout << "Exact solution: u(x) = sin(pi*x)\n";
    std::cout << "Newton status: " << result.nonlinear_solver.message << '\n';
    std::cout << "Newton iterations: " << result.nonlinear_solver.iterations << '\n';
    std::cout << "Residual norm: " << std::scientific
              << result.nonlinear_solver.residual_norm << '\n';

    if (!result.nonlinear_solver.converged) {
        std::cerr << "Nonlinear solve did not converge.\n";
        return 2;
    }

    double maximum_nodal_error = 0.0;
    std::cout << "\n       x          FEM u(x)       exact u(x)      abs. error\n";
    std::cout << std::fixed << std::setprecision(8);

    for (std::size_t node = 0; node < result.coordinates.size(); ++node) {
        const double x = result.coordinates[node];
        const double numerical = result.nonlinear_solver.solution[node];
        const double exact = exact_solution(x);
        const double error = std::abs(numerical - exact);
        maximum_nodal_error = std::max(maximum_nodal_error, error);

        std::cout << std::setw(10) << x << "  "
                  << std::setw(14) << numerical << "  "
                  << std::setw(14) << exact << "  "
                  << std::scientific << error << std::fixed << '\n';
    }

    std::cout << "\nMaximum nodal error: " << std::scientific
              << maximum_nodal_error << '\n';
    return 0;
}
