#include "test_framework.hpp"

#include "nova/fem/diffusion_2d.hpp"
#include "nova/fem/mesh2d.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

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

double manufactured_solution(const double x, const double y) {
    return x * (1.0 - x) * y * (1.0 - y);
}

double manufactured_source(const double x, const double y) {
    return 2.0 * (x * (1.0 - x) + y * (1.0 - y));
}

double solve_manufactured_problem(const std::size_t subdivisions) {
    auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, subdivisions, subdivisions
    );

    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(generated.mesh);
    problem.diffusivity = [](const double, const double) { return 1.0; };
    problem.source = manufactured_source;

    for (const std::size_t node : generated.boundary_nodes) {
        const auto& point = problem.mesh.nodes[node];
        problem.dirichlet_constraints.push_back({
            node,
            manufactured_solution(point.x, point.y)
        });
    }

    const auto result = nova::fem::solve_scalar_diffusion_2d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);

    double maximum_nodal_error = 0.0;
    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        const auto& point = problem.mesh.nodes[node];
        const double exact = manufactured_solution(point.x, point.y);
        const double error = std::abs(result.values[node] - exact);
        if (error > maximum_nodal_error) {
            maximum_nodal_error = error;
        }
    }
    return maximum_nodal_error;
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

NOVA_TEST(diffusion_2d_supports_variable_diffusivity_and_source) {
    auto problem = make_harmonic_problem();
    problem.diffusivity = [](const double x, const double y) {
        return 1.0 + x + y;
    };
    problem.source = [](const double, const double) { return -1.0; };
    problem.dirichlet_constraints = {
        {0, 0.0},
        {1, 1.0},
        {2, 1.0},
        {3, 0.0}
    };

    const auto result = nova::fem::solve_scalar_diffusion_2d(problem);

    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.values[4], 0.5, 1.0e-10);
}

NOVA_TEST(diffusion_2d_rejects_non_positive_diffusivity) {
    auto problem = make_harmonic_problem();
    problem.diffusivity = [](const double, const double) { return 0.0; };

    bool threw = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_2d(problem));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
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

NOVA_TEST(rectangle_mesh_generator_creates_expected_topology) {
    const auto generated = nova::fem::make_rectangle_mesh_2d(
        -1.0, 2.0, 0.0, 1.0, 4, 3
    );

    NOVA_REQUIRE(generated.mesh.nodes.size() == 20);
    NOVA_REQUIRE(generated.mesh.elements.size() == 24);
    NOVA_REQUIRE(generated.boundary_nodes.size() == 14);
    NOVA_REQUIRE(generated.boundary_edges.size() == 14);
    generated.mesh.validate();
}

NOVA_TEST(rectangle_mesh_generator_rejects_zero_subdivisions) {
    bool threw = false;
    try {
        static_cast<void>(nova::fem::make_rectangle_mesh_2d(
            0.0, 1.0, 0.0, 1.0, 0, 4
        ));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(diffusion_2d_solves_mixed_dirichlet_neumann_problem) {
    auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 4, 4
    );

    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(generated.mesh);
    problem.diffusivity = [](const double, const double) { return 1.0; };
    problem.source = [](const double, const double) { return 0.0; };

    for (const std::size_t node : generated.boundary_nodes) {
        const auto& point = problem.mesh.nodes[node];
        if (point.x == 0.0) {
            problem.dirichlet_constraints.push_back({node, 0.0});
        }
    }

    for (const auto& edge : generated.boundary_edges) {
        const auto& p0 = problem.mesh.nodes[edge[0]];
        const auto& p1 = problem.mesh.nodes[edge[1]];
        if (p0.x == 0.0 && p1.x == 0.0) {
            continue;
        }

        const double flux = (p0.x == 1.0 && p1.x == 1.0) ? -1.0 : 0.0;
        problem.neumann_conditions.push_back({
            edge,
            [flux](const double, const double) { return flux; }
        });
    }

    const auto result = nova::fem::solve_scalar_diffusion_2d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);

    for (std::size_t node = 0; node < problem.mesh.nodes.size(); ++node) {
        NOVA_REQUIRE_NEAR(
            result.values[node],
            problem.mesh.nodes[node].x,
            1.0e-9
        );
    }
}

NOVA_TEST(diffusion_2d_rejects_neumann_condition_on_interior_edge) {
    auto problem = make_harmonic_problem();
    problem.neumann_conditions.push_back({
        {0, 4},
        [](const double, const double) { return 0.0; }
    });

    bool threw = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_2d(problem));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(diffusion_2d_shows_second_order_nodal_convergence) {
    const double coarse_error = solve_manufactured_problem(4);
    const double medium_error = solve_manufactured_problem(8);
    const double fine_error = solve_manufactured_problem(16);

    NOVA_REQUIRE(coarse_error > medium_error);
    NOVA_REQUIRE(medium_error > fine_error);
    NOVA_REQUIRE(coarse_error / medium_error > 3.5);
    NOVA_REQUIRE(medium_error / fine_error > 3.5);
}
