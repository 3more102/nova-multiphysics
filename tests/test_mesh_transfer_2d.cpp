#include "test_framework.hpp"

#include "nova/fem/mesh_transfer_2d.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using nova::fem::H1FunctionSpace2D;
using nova::fem::H1TriangleOrder;
using nova::fem::TransientHeat2DProblem;
using nova::fem::TriangleMesh2D;

void check_exact_uniform_transfer(
    const H1TriangleOrder order,
    const std::function<double(double, double)>& field
) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        -1.0, 2.0, -0.5, 1.5, 2U, 3U
    );
    const H1FunctionSpace2D coarse_space(rectangle.mesh, order);
    const auto coarse_values = coarse_space.interpolate(field);
    const auto refined = nova::fem::refine_triangles_uniformly_2d(
        rectangle.mesh, rectangle.boundary_edges
    );

    const auto fine_values = nova::fem::prolongate_h1_uniform_2d(
        rectangle.mesh, refined, order, coarse_values
    );
    const H1FunctionSpace2D fine_space(refined.mesh, order);
    const auto expected = fine_space.interpolate(field);

    NOVA_REQUIRE(fine_values.size() == fine_space.dimension());
    NOVA_REQUIRE(expected.size() == fine_values.size());
    for (std::size_t i = 0; i < fine_values.size(); ++i) {
        NOVA_REQUIRE_NEAR(fine_values[i], expected[i], 2.0e-12);
    }
}

} // namespace

NOVA_TEST(mesh_transfer_preserves_affine_p1_field_exactly) {
    check_exact_uniform_transfer(
        H1TriangleOrder::P1,
        [](const double x, const double y) {
            return 1.25 + 2.0 * x - 3.0 * y;
        }
    );
}

NOVA_TEST(mesh_transfer_preserves_quadratic_p2_field_exactly) {
    check_exact_uniform_transfer(
        H1TriangleOrder::P2,
        [](const double x, const double y) {
            return 1.0 + 2.0 * x - 3.0 * y +
                0.5 * x * x + x * y + 2.0 * y * y;
        }
    );
}

NOVA_TEST(mesh_transfer_rejects_invalid_state_and_parent_map) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    const auto refined = nova::fem::refine_triangles_uniformly_2d(
        rectangle.mesh
    );
    const H1FunctionSpace2D coarse_space(
        rectangle.mesh, H1TriangleOrder::P2
    );
    auto values = coarse_space.interpolate(
        [](const double x, const double y) { return x + y; }
    );

    bool rejected_size = false;
    try {
        static_cast<void>(nova::fem::prolongate_h1_uniform_2d(
            rectangle.mesh, refined, H1TriangleOrder::P2,
            std::vector<double>(1U, 0.0)
        ));
    } catch (const std::invalid_argument&) {
        rejected_size = true;
    }
    NOVA_REQUIRE(rejected_size);

    values[0] = std::numeric_limits<double>::quiet_NaN();
    bool rejected_nonfinite = false;
    try {
        static_cast<void>(nova::fem::prolongate_h1_uniform_2d(
            rectangle.mesh, refined, H1TriangleOrder::P2, values
        ));
    } catch (const std::invalid_argument&) {
        rejected_nonfinite = true;
    }
    NOVA_REQUIRE(rejected_nonfinite);

    values[0] = 0.0;
    auto invalid_map = refined;
    invalid_map.parent_elements[0] = rectangle.mesh.elements.size();
    bool rejected_parent = false;
    try {
        static_cast<void>(nova::fem::prolongate_h1_uniform_2d(
            rectangle.mesh, invalid_map, H1TriangleOrder::P2, values
        ));
    } catch (const std::invalid_argument&) {
        rejected_parent = true;
    }
    NOVA_REQUIRE(rejected_parent);
}

NOVA_TEST(transient_heat_uniform_refinement_transfers_p2_state_and_boundaries) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );

    TransientHeat2DProblem problem;
    problem.mesh = rectangle.mesh;
    problem.volumetric_capacity = [](const double, const double) { return 1.0; };
    problem.conductivity = [](const double, const double) { return 1.0; };
    problem.source = [](const double, const double, const double) { return 0.0; };
    problem.initial_temperature = [](const double, const double) { return 0.0; };

    const H1FunctionSpace2D coarse_space(problem.mesh, H1TriangleOrder::P2);
    problem.initial_dof_values = coarse_space.interpolate(
        [](const double x, const double y) {
            return 0.5 + x - 2.0 * y + x * x + x * y;
        }
    );

    problem.dirichlet_edges.push_back({
        rectangle.boundary_edges[0],
        [](const double, const double, const double) { return 0.0; }
    });
    problem.neumann_edges.push_back({
        rectangle.boundary_edges[1],
        [](const double, const double, const double) { return 0.0; }
    });
    problem.robin_edges.push_back({
        rectangle.boundary_edges[2],
        [](const double, const double, const double) { return 0.0; },
        [](const double, const double, const double) { return 0.0; }
    });

    const auto refined =
        nova::fem::refine_transient_heat_problem_uniformly_2d(
            problem, H1TriangleOrder::P2
        );

    NOVA_REQUIRE(refined.problem.mesh.elements.size() ==
        4U * problem.mesh.elements.size());
    NOVA_REQUIRE(refined.problem.dirichlet_edges.size() == 2U);
    NOVA_REQUIRE(refined.problem.neumann_edges.size() == 2U);
    NOVA_REQUIRE(refined.problem.robin_edges.size() == 2U);

    const H1FunctionSpace2D fine_space(
        refined.problem.mesh, H1TriangleOrder::P2
    );
    NOVA_REQUIRE(refined.problem.initial_dof_values.size() ==
        fine_space.dimension());

    for (const double value : refined.problem.initial_dof_values) {
        NOVA_REQUIRE(std::isfinite(value));
    }

    nova::fem::TransientHeat2DOptions options;
    options.order = H1TriangleOrder::P2;
    options.time_step = 0.1;
    options.steps = 1U;
    const auto solution = nova::fem::solve_transient_heat_2d(
        refined.problem, options
    );
    NOVA_REQUIRE(solution.temperature.size() == fine_space.dimension());
    for (const double value : solution.temperature) {
        NOVA_REQUIRE(std::isfinite(value));
    }
}

NOVA_TEST(transient_heat_uniform_refinement_rejects_duplicate_boundary_assignment) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    TransientHeat2DProblem problem;
    problem.mesh = rectangle.mesh;
    problem.dirichlet_edges.push_back({
        rectangle.boundary_edges[0],
        [](const double, const double, const double) { return 0.0; }
    });
    problem.neumann_edges.push_back({
        {rectangle.boundary_edges[0][1], rectangle.boundary_edges[0][0]},
        [](const double, const double, const double) { return 0.0; }
    });

    bool rejected = false;
    try {
        static_cast<void>(
            nova::fem::refine_transient_heat_problem_uniformly_2d(
                problem, H1TriangleOrder::P1
            )
        );
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}


NOVA_TEST(mesh_transfer_handles_large_coordinate_offsets) {
    const TriangleMesh2D mesh{{
        {1000000.1, 1000000.2},
        {1000000.4, 1000000.3},
        {1000000.2, 1000000.6}
    }, {{0U, 1U, 2U}}};
    const auto refined = nova::fem::refine_triangles_uniformly_2d(mesh);
    const H1FunctionSpace2D coarse_space(mesh, H1TriangleOrder::P2);
    const auto coarse_values = coarse_space.interpolate(
        [](const double, const double) { return 1.0; }
    );
    const auto fine_values = nova::fem::prolongate_h1_uniform_2d(
        mesh, refined, H1TriangleOrder::P2, coarse_values
    );
    for (const double value : fine_values) {
        NOVA_REQUIRE_NEAR(value, 1.0, 2.0e-13);
    }
}
