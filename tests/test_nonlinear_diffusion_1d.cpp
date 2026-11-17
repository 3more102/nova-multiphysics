#include "test_framework.hpp"

#include "nova/fem/nonlinear_diffusion_1d.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

constexpr double pi = 3.14159265358979323846;

double exact_solution(const double x) {
    return std::sin(pi * x);
}

double manufactured_source(const double x) {
    const double u = exact_solution(x);
    const double derivative = pi * std::cos(pi * x);
    const double second_derivative = -pi * pi * u;

    return -(
        2.0 * u * derivative * derivative +
        (1.0 + u * u) * second_derivative
    );
}

double maximum_nodal_error(const std::size_t elements) {
    nova::fem::NonlinearDiffusion1DProblem problem;
    problem.length = 1.0;
    problem.elements = elements;
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

    NOVA_REQUIRE(result.nonlinear_solver.converged);
    NOVA_REQUIRE(result.nonlinear_solver.iterations > 0);

    double error = 0.0;
    for (std::size_t node = 0; node < result.coordinates.size(); ++node) {
        const double exact = exact_solution(result.coordinates[node]);
        const double difference =
            std::abs(result.nonlinear_solver.solution[node] - exact);
        if (difference > error) {
            error = difference;
        }
    }
    return error;
}

} // namespace

NOVA_TEST(nonlinear_diffusion_1d_converges_to_manufactured_solution) {
    const double error_8 = maximum_nodal_error(8);
    const double error_16 = maximum_nodal_error(16);
    const double error_32 = maximum_nodal_error(32);

    NOVA_REQUIRE(error_8 > error_16);
    NOVA_REQUIRE(error_16 > error_32);
    NOVA_REQUIRE(error_8 < 1.0e-5);
    NOVA_REQUIRE(error_32 < 1.0e-7);
}

NOVA_TEST(nonlinear_diffusion_1d_supports_nonzero_dirichlet_values) {
    nova::fem::NonlinearDiffusion1DProblem problem;
    problem.length = 1.0;
    problem.elements = 8;
    problem.source = [](const double) { return 0.0; };
    problem.left_value = 1.0;
    problem.right_value = 2.0;

    const auto result = nova::fem::solve_nonlinear_diffusion_1d(problem);

    NOVA_REQUIRE(result.nonlinear_solver.converged);
    NOVA_REQUIRE_NEAR(result.nonlinear_solver.solution.front(), 1.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.nonlinear_solver.solution.back(), 2.0, 1.0e-12);
}

NOVA_TEST(nonlinear_diffusion_1d_rejects_invalid_initial_guess_size) {
    nova::fem::NonlinearDiffusion1DProblem problem;
    problem.elements = 4;
    problem.source = [](const double) { return 0.0; };
    const std::vector<double> initial_guess{0.0, 0.0};

    bool threw = false;
    try {
        static_cast<void>(nova::fem::solve_nonlinear_diffusion_1d(
            problem,
            initial_guess
        ));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
