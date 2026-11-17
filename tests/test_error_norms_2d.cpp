#include "test_framework.hpp"

#include "nova/fem/diffusion_2d.hpp"
#include "nova/fem/error_norms_2d.hpp"
#include "nova/fem/mesh2d.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

double exact_value(const double x, const double y) {
    return x * (1.0 - x) * y * (1.0 - y);
}

double exact_derivative_x(const double x, const double y) {
    return (1.0 - 2.0 * x) * y * (1.0 - y);
}

double exact_derivative_y(const double x, const double y) {
    return x * (1.0 - x) * (1.0 - 2.0 * y);
}

nova::fem::ErrorNorms2D solve_and_measure(const std::size_t subdivisions) {
    auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, subdivisions, subdivisions
    );

    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = generated.mesh;
    problem.diffusivity = [](const double, const double) { return 1.0; };
    problem.source = [](const double x, const double y) {
        return 2.0 * (x * (1.0 - x) + y * (1.0 - y));
    };

    for (const std::size_t node : generated.boundary_nodes) {
        const auto& point = problem.mesh.nodes[node];
        problem.dirichlet_constraints.push_back({
            node,
            exact_value(point.x, point.y)
        });
    }

    const auto solution = nova::fem::solve_scalar_diffusion_2d(problem);
    NOVA_REQUIRE(solution.linear_solver.converged);

    const nova::fem::ExactScalarField2D exact{
        exact_value,
        exact_derivative_x,
        exact_derivative_y
    };

    return nova::fem::compute_error_norms_2d(
        problem.mesh,
        solution.values,
        exact
    );
}

} // namespace

NOVA_TEST(error_norms_are_zero_for_exact_linear_p1_field) {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
        {0.5, 0.5}
    };
    mesh.elements = {
        {0, 1, 4},
        {1, 2, 4},
        {2, 3, 4},
        {3, 0, 4}
    };

    std::vector<double> values;
    values.reserve(mesh.nodes.size());
    for (const auto& point : mesh.nodes) {
        values.push_back(1.0 + 2.0 * point.x - 3.0 * point.y);
    }

    const nova::fem::ExactScalarField2D exact{
        [](const double x, const double y) { return 1.0 + 2.0 * x - 3.0 * y; },
        [](const double, const double) { return 2.0; },
        [](const double, const double) { return -3.0; }
    };

    const auto errors = nova::fem::compute_error_norms_2d(mesh, values, exact);
    NOVA_REQUIRE(errors.l2 < 1.0e-13);
    NOVA_REQUIRE(errors.h1_seminorm < 1.0e-13);
    NOVA_REQUIRE(errors.h1 < 1.0e-13);
}

NOVA_TEST(error_norms_are_orientation_independent) {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
        {0.5, 0.5}
    };
    mesh.elements = {
        {4, 1, 0},
        {1, 2, 4},
        {2, 3, 4},
        {3, 0, 4}
    };

    std::vector<double> values;
    values.reserve(mesh.nodes.size());
    for (const auto& point : mesh.nodes) {
        values.push_back(1.0 + 2.0 * point.x - 3.0 * point.y);
    }

    const nova::fem::ExactScalarField2D exact{
        [](const double x, const double y) { return 1.0 + 2.0 * x - 3.0 * y; },
        [](const double, const double) { return 2.0; },
        [](const double, const double) { return -3.0; }
    };

    const auto errors = nova::fem::compute_error_norms_2d(mesh, values, exact);
    NOVA_REQUIRE(errors.l2 < 1.0e-13);
    NOVA_REQUIRE(errors.h1_seminorm < 1.0e-13);
}

NOVA_TEST(error_norms_converge_for_manufactured_quartic_solution) {
    const auto coarse = solve_and_measure(4);
    const auto medium = solve_and_measure(8);
    const auto fine = solve_and_measure(16);

    NOVA_REQUIRE(coarse.l2 > medium.l2);
    NOVA_REQUIRE(medium.l2 > fine.l2);
    NOVA_REQUIRE(coarse.l2 / medium.l2 > 3.5);
    NOVA_REQUIRE(medium.l2 / fine.l2 > 3.5);

    NOVA_REQUIRE(coarse.h1_seminorm > medium.h1_seminorm);
    NOVA_REQUIRE(medium.h1_seminorm > fine.h1_seminorm);
    NOVA_REQUIRE(coarse.h1_seminorm / medium.h1_seminorm > 1.7);
    NOVA_REQUIRE(medium.h1_seminorm / fine.h1_seminorm > 1.7);

    NOVA_REQUIRE(coarse.h1 > medium.h1);
    NOVA_REQUIRE(medium.h1 > fine.h1);
}

NOVA_TEST(error_norms_reject_field_size_mismatch) {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    mesh.elements = {{0, 1, 2}};

    const std::vector<double> values{0.0, 0.0};
    const nova::fem::ExactScalarField2D exact{
        [](const double, const double) { return 0.0; },
        [](const double, const double) { return 0.0; },
        [](const double, const double) { return 0.0; }
    };

    bool threw = false;
    try {
        static_cast<void>(nova::fem::compute_error_norms_2d(mesh, values, exact));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
