#include "nova/fem/poisson_1d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    nova::fem::Poisson1DProblem problem;
    problem.length = 1.0;
    problem.elements = 16;
    problem.conductivity = 1.0;
    problem.source = [](const double) { return 2.0; };
    problem.left_value = 0.0;
    problem.right_value = 0.0;

    const nova::fem::Poisson1DResult result = nova::fem::solve_poisson_1d(problem);

    std::cout << "NOVA Multiphysics - 1D Poisson benchmark\n";
    std::cout << "Equation: -u'' = 2 on [0, 1], with u(0) = u(1) = 0\n";
    std::cout << "Analytical solution: u(x) = x(1-x)\n";
    std::cout << "Solver: " << result.linear_solver.message << '\n';
    std::cout << "Iterations: " << result.linear_solver.iterations << '\n';
    std::cout << "Residual norm: " << std::scientific
              << result.linear_solver.residual_norm << "\n\n";

    if (!result.linear_solver.converged) {
        std::cerr << "The numerical solve did not converge.\n";
        return 2;
    }

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "       x          FEM u(x)       exact u(x)      abs. error\n";

    double maximum_nodal_error = 0.0;
    for (std::size_t index = 0; index < result.coordinates.size(); ++index) {
        const double x = result.coordinates[index];
        const double exact = x * (1.0 - x);
        const double error = std::abs(result.values[index] - exact);
        maximum_nodal_error = std::max(maximum_nodal_error, error);

        std::cout << std::setw(10) << x << "  "
                  << std::setw(14) << result.values[index] << "  "
                  << std::setw(14) << exact << "  "
                  << std::scientific << error << std::fixed << '\n';
    }

    std::cout << "\nMaximum nodal error: " << std::scientific
              << maximum_nodal_error << '\n';

    return 0;
}
