#include "test_framework.hpp"

#include "nova/fem/mesh_adaptivity_2d.hpp"
#include "nova/fem/mesh_transfer_2d.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using nova::fem::BoundaryEdge2D;
using nova::fem::H1FunctionSpace2D;
using nova::fem::H1TriangleOrder;
using nova::fem::TriangleMesh2D;

[[nodiscard]] std::pair<std::size_t, std::size_t> edge_key(
    const std::size_t a, const std::size_t b
) {
    return std::minmax(a, b);
}

[[nodiscard]] std::map<std::pair<std::size_t, std::size_t>, std::size_t>
edge_counts(const TriangleMesh2D& mesh) {
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> counts;
    for (const auto& triangle : mesh.elements) {
        ++counts[edge_key(triangle[0], triangle[1])];
        ++counts[edge_key(triangle[1], triangle[2])];
        ++counts[edge_key(triangle[2], triangle[0])];
    }
    return counts;
}

void check_affine_recovery(const H1TriangleOrder order) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        -1.0, 2.0, -0.5, 1.5, 3U, 2U
    );
    const H1FunctionSpace2D space(rectangle.mesh, order);
    const auto values = space.interpolate([](const double x, const double y) {
        return 1.25 + 2.0 * x - 3.0 * y;
    });
    const auto indicators =
        nova::fem::estimate_h1_gradient_recovery_indicators_2d(
            rectangle.mesh, order, values
        );
    NOVA_REQUIRE(indicators.size() == rectangle.mesh.elements.size());
    for (const double indicator : indicators) {
        NOVA_REQUIRE(indicator < 1.0e-12);
    }
}

} // namespace

NOVA_TEST(mesh_adaptivity_gradient_recovery_is_zero_for_affine_p1_and_p2) {
    check_affine_recovery(H1TriangleOrder::P1);
    check_affine_recovery(H1TriangleOrder::P2);
}

NOVA_TEST(mesh_adaptivity_gradient_recovery_detects_curved_fields) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 4U, 4U
    );
    for (const H1TriangleOrder order : {
             H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        const H1FunctionSpace2D space(rectangle.mesh, order);
        const auto values = space.interpolate([](const double x, const double y) {
            return std::sin(3.0 * x) * std::cos(2.0 * y);
        });
        const auto indicators =
            nova::fem::estimate_h1_gradient_recovery_indicators_2d(
                rectangle.mesh, order, values
            );
        double squared_sum = 0.0;
        for (const double indicator : indicators) {
            NOVA_REQUIRE(std::isfinite(indicator));
            NOVA_REQUIRE(indicator >= 0.0);
            squared_sum += indicator * indicator;
        }
        NOVA_REQUIRE(squared_sum > 1.0e-12);
    }
}

NOVA_TEST(mesh_adaptivity_dorfler_bulk_marking_is_minimal_and_deterministic) {
    const std::vector<double> indicators{1.0, 2.0, 3.0, 4.0};
    const auto half = nova::fem::mark_elements_dorfler_2d(indicators, 0.5);
    NOVA_REQUIRE(half.size() == 1U);
    NOVA_REQUIRE(half[0] == 3U);

    const auto eighty = nova::fem::mark_elements_dorfler_2d(indicators, 0.8);
    NOVA_REQUIRE(eighty.size() == 2U);
    NOVA_REQUIRE(eighty[0] == 3U);
    NOVA_REQUIRE(eighty[1] == 2U);

    const auto zeros = nova::fem::mark_elements_dorfler_2d(
        std::vector<double>{0.0, 0.0}, 0.5
    );
    NOVA_REQUIRE(zeros.empty());

    bool rejected_theta = false;
    try {
        static_cast<void>(nova::fem::mark_elements_dorfler_2d(indicators, 0.0));
    } catch (const std::invalid_argument&) {
        rejected_theta = true;
    }
    NOVA_REQUIRE(rejected_theta);

    bool rejected_negative = false;
    try {
        static_cast<void>(nova::fem::mark_elements_dorfler_2d(
            std::vector<double>{1.0, -1.0}, 0.5
        ));
    } catch (const std::invalid_argument&) {
        rejected_negative = true;
    }
    NOVA_REQUIRE(rejected_negative);
}

NOVA_TEST(mesh_adaptivity_local_refinement_closes_shared_edges_without_hanging_nodes) {
    TriangleMesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}
    };
    mesh.elements = {{0U, 1U, 2U}, {0U, 2U, 3U}};
    const std::vector<BoundaryEdge2D> boundary{
        {0U, 1U}, {1U, 2U}, {2U, 3U}, {3U, 0U}
    };
    const auto old_quality = nova::fem::evaluate_mesh_quality_2d(mesh);
    const auto refined =
        nova::fem::refine_triangles_by_marked_elements_2d(
            mesh, boundary, std::vector<std::size_t>{0U}
        );

    NOVA_REQUIRE(refined.midpoint_parent_edges.size() == 1U);
    NOVA_REQUIRE(refined.midpoint_parent_edges[0] == BoundaryEdge2D{0U, 2U});
    NOVA_REQUIRE(refined.mesh.nodes.size() == mesh.nodes.size() + 1U);
    NOVA_REQUIRE(refined.mesh.elements.size() == 4U);
    NOVA_REQUIRE(refined.parent_elements.size() == 4U);
    NOVA_REQUIRE(refined.boundary_edges.size() == boundary.size());

    const auto counts = edge_counts(refined.mesh);
    for (const auto& [edge, count] : counts) {
        static_cast<void>(edge);
        NOVA_REQUIRE(count == 1U || count == 2U);
    }
    const auto new_quality = nova::fem::evaluate_mesh_quality_2d(refined.mesh);
    NOVA_REQUIRE_NEAR(old_quality.total_area, new_quality.total_area, 1.0e-13);
    for (const auto& quality : new_quality.elements) {
        NOVA_REQUIRE(quality.signed_area > 0.0);
    }
}

NOVA_TEST(mesh_adaptivity_local_boundary_refinement_preserves_direction) {
    TriangleMesh2D mesh;
    mesh.nodes = {{0.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}};
    mesh.elements = {{0U, 1U, 2U}};
    const std::vector<BoundaryEdge2D> boundary{
        {0U, 1U}, {1U, 2U}, {2U, 0U}
    };
    const auto refined =
        nova::fem::refine_triangles_by_marked_elements_2d(
            mesh, boundary, std::vector<std::size_t>{0U}
        );
    NOVA_REQUIRE(refined.midpoint_parent_edges.size() == 1U);
    NOVA_REQUIRE(refined.midpoint_parent_edges[0] == BoundaryEdge2D{1U, 2U});
    NOVA_REQUIRE(refined.boundary_edges.size() == 4U);
    NOVA_REQUIRE(refined.boundary_edges[0] == BoundaryEdge2D{0U, 1U});
    NOVA_REQUIRE(refined.boundary_edges[1][0] == 1U);
    NOVA_REQUIRE(refined.boundary_edges[1][1] == 3U);
    NOVA_REQUIRE(refined.boundary_edges[2][0] == 3U);
    NOVA_REQUIRE(refined.boundary_edges[2][1] == 2U);
    NOVA_REQUIRE(refined.boundary_edges[3] == BoundaryEdge2D{2U, 0U});
    const auto quality = nova::fem::evaluate_mesh_quality_2d(refined.mesh);
    NOVA_REQUIRE_NEAR(quality.total_area, 1.0, 1.0e-13);
}

NOVA_TEST(mesh_adaptivity_rejects_duplicate_and_out_of_range_marks) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    bool duplicate = false;
    try {
        static_cast<void>(nova::fem::refine_triangles_by_marked_elements_2d(
            rectangle.mesh, {}, std::vector<std::size_t>{0U, 0U}
        ));
    } catch (const std::invalid_argument&) {
        duplicate = true;
    }
    NOVA_REQUIRE(duplicate);

    bool out_of_range = false;
    try {
        static_cast<void>(nova::fem::refine_triangles_by_marked_elements_2d(
            rectangle.mesh, {}, std::vector<std::size_t>{99U}
        ));
    } catch (const std::out_of_range&) {
        out_of_range = true;
    }
    NOVA_REQUIRE(out_of_range);
}

NOVA_TEST(transient_heat_adaptive_remeshing_transfers_p2_and_boundary_callbacks) {
    TriangleMesh2D mesh;
    mesh.nodes = {{0.0, 0.0}, {2.0, 0.0}, {0.0, 1.0}};
    mesh.elements = {{0U, 1U, 2U}};

    nova::fem::TransientHeat2DProblem problem;
    problem.mesh = mesh;
    problem.volumetric_capacity =
        [](const double, const double) { return 1.0; };
    problem.conductivity =
        [](const double, const double) { return 1.0; };
    problem.source =
        [](const double, const double, const double) { return 0.0; };
    problem.initial_temperature =
        [](const double x, const double y) { return 0.0; };

    const H1FunctionSpace2D coarse_space(mesh, H1TriangleOrder::P2);
    const auto exact_field = [](const double x, const double y) {
        return 0.5 + x - 2.0 * y + x * x + x * y;
    };
    problem.initial_dof_values = coarse_space.interpolate(exact_field);
    problem.dirichlet_edges.push_back({
        {1U, 2U},
        [](const double x, const double y, const double t) {
            return 0.5 + x - 2.0 * y + x * x + x * y + t;
        }
    });
    problem.neumann_edges.push_back({
        {2U, 0U},
        [](const double, const double, const double) { return 0.0; }
    });
    problem.robin_edges.push_back({
        {0U, 1U},
        [](const double, const double, const double) { return 1.0; },
        [](const double, const double, const double) { return 0.0; }
    });

    const auto refined =
        nova::fem::refine_transient_heat_problem_adaptively_2d(
            problem, H1TriangleOrder::P2, std::vector<std::size_t>{0U}
        );
    NOVA_REQUIRE(refined.problem.mesh.elements.size() == 2U);
    NOVA_REQUIRE(refined.problem.dirichlet_edges.size() == 2U);
    NOVA_REQUIRE(refined.problem.neumann_edges.size() == 1U);
    NOVA_REQUIRE(refined.problem.robin_edges.size() == 1U);

    const H1FunctionSpace2D fine_space(
        refined.problem.mesh, H1TriangleOrder::P2
    );
    const auto expected = fine_space.interpolate(exact_field);
    NOVA_REQUIRE(refined.problem.initial_dof_values.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        NOVA_REQUIRE_NEAR(
            refined.problem.initial_dof_values[i], expected[i], 2.0e-12
        );
    }

    nova::fem::TransientHeat2DOptions options;
    options.order = H1TriangleOrder::P2;
    options.time_step = 0.01;
    options.steps = 1U;
    const auto solution = nova::fem::solve_transient_heat_2d(
        refined.problem, options
    );
    NOVA_REQUIRE(solution.temperature.size() == fine_space.dimension());
    for (const double value : solution.temperature) {
        NOVA_REQUIRE(std::isfinite(value));
    }
}
