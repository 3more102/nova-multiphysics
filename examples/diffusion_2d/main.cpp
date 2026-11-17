#include "nova/fem/diffusion_2d.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
        {0.5, 0.5}
    };
    problem.mesh.elements = {
        {0, 1, 4},
        {1, 2, 4},
        {2, 3, 4},
        {3, 0, 4}
    };

    // Exact harmonic solution: u(x,y) = 1 + 2x - 3y.
    problem.diffusivity = [](const double, const double) { return 1.0; };
    problem.source = [](const double, const double) { return 0.0; };
    problem.dirichlet_constraints = {
        {0, 1.0},
        {1, 3.0},
        {2, 0.0},
        {3, -2.0}
    };

    const auto result = nova::fem::solve_scalar_diffusion_2d(problem);

    std::cout << "NOVA - 2D P1 triangular diffusion benchmark\n";
    std::cout << "Exact solution: u(x,y) = 1 + 2x - 3y\n";
    std::cout << "Solver: " << result.linear_solver.message << '\n';
    std::cout << "Iterations: " << result.linear_solver.iterations << '\n';
    std::cout << "Residual norm: " << std::scientific
              << result.linear_solver.residual_norm << '\n';

    if (!result.linear_solver.converged) {
        std::cerr << "The numerical solve did not converge.\n";
        return 2;
    }

    std::cout << std::fixed << std::setprecision(10);
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        const auto& point = problem.mesh.nodes[node];
        const double exact = 1.0 + 2.0 * point.x - 3.0 * point.y;
        std::cout << "node " << node
                  << "  FEM=" << result.values[node]
                  << "  exact=" << exact
                  << "  error=" << std::scientific
                  << std::abs(result.values[node] - exact)
                  << std::fixed << '\n';
    }

    return 0;
}
