#include "test_framework.hpp"

#include "nova/fem/poisson_1d.hpp"

#include <stdexcept>

NOVA_TEST(poisson_1d_matches_quadratic_analytical_solution_at_nodes) {
    nova::fem::Poisson1DProblem problem;
    problem.length = 1.0;
    problem.elements = 12;
    problem.conductivity = 1.0;
    problem.source = [](const double) { return 2.0; };
    problem.left_value = 0.0;
    problem.right_value = 0.0;

    const auto result = nova::fem::solve_poisson_1d(problem);

    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE(result.coordinates.size() == 13);
    NOVA_REQUIRE(result.values.size() == 13);

    for (std::size_t index = 0; index < result.coordinates.size(); ++index) {
        const double x = result.coordinates[index];
        const double exact = x * (1.0 - x);
        NOVA_REQUIRE_NEAR(result.values[index], exact, 1.0e-10);
    }
}

NOVA_TEST(poisson_1d_supports_nonzero_dirichlet_values) {
    nova::fem::Poisson1DProblem problem;
    problem.length = 2.0;
    problem.elements = 8;
    problem.conductivity = 3.0;
    problem.source = [](const double) { return 0.0; };
    problem.left_value = 2.0;
    problem.right_value = 5.0;

    const auto result = nova::fem::solve_poisson_1d(problem);

    NOVA_REQUIRE(result.linear_solver.converged);
    for (std::size_t index = 0; index < result.coordinates.size(); ++index) {
        const double x = result.coordinates[index];
        const double exact = 2.0 + 1.5 * x;
        NOVA_REQUIRE_NEAR(result.values[index], exact, 1.0e-10);
    }
}

NOVA_TEST(poisson_1d_rejects_zero_elements) {
    nova::fem::Poisson1DProblem problem;
    problem.elements = 0;
    problem.source = [](const double) { return 1.0; };

    bool threw = false;
    try {
        static_cast<void>(nova::fem::solve_poisson_1d(problem));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
