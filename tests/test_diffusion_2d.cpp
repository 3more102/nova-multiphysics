#include "test_framework.hpp"

#include "nova/fem/diffusion_2d.hpp"

#include <stdexcept>

namespace {

nova::fem::ScalarDiffusion2DProblem make_harmonic_problem() {
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
    problem.diffusivity = [](const double, const double) { return 2.0; };
    problem.source = [](const double, const double) { return 0.0; };
    problem.dirichlet_constraints = {
        {0, 1.0},
        {1, 3.0},
        {2, 0.0},
        {3, -2.0}
    };
    return problem;
}

} // namespace

NOVA_TEST(diffusion_2d_reproduces_linear_harmonic_field) {
    const auto problem = make_harmonic_problem();
    const auto result = nova::fem::solve_scalar_diffusion_2d(problem);

    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE(result.values.size() == problem.mesh.nodes.size());

    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        const auto& point = problem.mesh.nodes[node];
        const double exact = 1.0 + 2.0 * point.x - 3.0 * point.y;
        NOVA_REQUIRE_NEAR(result.values[node], exact, 1.0e-10);
    }
}

NOVA_TEST(diffusion_2d_accepts_reversed_triangle_orientation) {
    auto problem = make_harmonic_problem();
    problem.mesh.elements[0] = {4, 1, 0};

    const auto result = nova::fem::solve_scalar_diffusion_2d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.values[4], 0.5, 1.0e-10);
}

NOVA_TEST(diffusion_2d_rejects_degenerate_triangle) {
    auto problem = make_harmonic_problem();
    problem.mesh.nodes[4] = {0.5, 0.0};

    bool threw = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_2d(problem));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(diffusion_2d_rejects_duplicate_dirichlet_constraint) {
    auto problem = make_harmonic_problem();
    problem.dirichlet_constraints.push_back({0, 1.0});

    bool threw = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_2d(problem));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
