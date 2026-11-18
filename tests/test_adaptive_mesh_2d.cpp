#include "test_framework.hpp"

#include "nova/fem/mesh_transfer_2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace {

using nova::fem::BoundaryEdge2D;
using nova::fem::H1FunctionSpace2D;
using nova::fem::H1TriangleOrder;
using nova::fem::TriangleMesh2D;

using EdgeKey = std::pair<std::size_t, std::size_t>;

[[nodiscard]] EdgeKey key(const std::size_t a, const std::size_t b) {
    return std::minmax(a, b);
}

[[nodiscard]] std::map<EdgeKey, std::size_t> edge_counts(
    const TriangleMesh2D& mesh
) {
    std::map<EdgeKey, std::size_t> counts;
    for (const auto& triangle : mesh.elements) {
        ++counts[key(triangle[0], triangle[1])];
        ++counts[key(triangle[1], triangle[2])];
        ++counts[key(triangle[2], triangle[0])];
    }
    return counts;
}

} // namespace

NOVA_TEST(adaptive_refinement_red_green_closure_is_conforming_and_conservative) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    const std::array<std::size_t, 1> marked{{0U}};
    const auto refined = nova::fem::refine_triangles_adaptively_2d(
        rectangle.mesh, marked, rectangle.boundary_edges
    );

    NOVA_REQUIRE(refined.mesh.elements.size() == 6U);
    NOVA_REQUIRE(refined.parent_elements.size() == 6U);
    NOVA_REQUIRE(refined.coarse_node_count == rectangle.mesh.nodes.size());
    NOVA_REQUIRE(refined.midpoint_parent_edges.size() == 3U);

    std::size_t children_of_zero = 0U;
    std::size_t children_of_one = 0U;
    for (const std::size_t parent : refined.parent_elements) {
        if (parent == 0U) {
            ++children_of_zero;
        } else if (parent == 1U) {
            ++children_of_one;
        }
    }
    NOVA_REQUIRE(children_of_zero == 4U);
    NOVA_REQUIRE(children_of_one == 2U);

    const auto coarse_quality = nova::fem::evaluate_mesh_quality_2d(rectangle.mesh);
    const auto fine_quality = nova::fem::evaluate_mesh_quality_2d(refined.mesh);
    NOVA_REQUIRE_NEAR(
        fine_quality.total_area, coarse_quality.total_area, 1.0e-13
    );

    const auto counts = edge_counts(refined.mesh);
    for (const auto& [edge, incidence] : counts) {
        static_cast<void>(edge);
        NOVA_REQUIRE(incidence == 1U || incidence == 2U);
    }
    NOVA_REQUIRE(refined.boundary_edges.size() == 6U);
    for (const BoundaryEdge2D& edge : refined.boundary_edges) {
        NOVA_REQUIRE(counts.at(key(edge[0], edge[1])) == 1U);
    }
}

NOVA_TEST(adaptive_refinement_two_and_three_edge_green_transitions_are_conforming) {
    const double root_three = std::sqrt(3.0);
    TriangleMesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0}, {2.0, 0.0}, {1.0, root_three},
        {1.0, -0.2},
        {1.5 + 0.1 * root_three, root_three / 2.0 + 0.1},
        {0.5 - 0.1 * root_three, root_three / 2.0 + 0.1}
    };
    mesh.elements = {
        {0U, 1U, 2U},
        {1U, 0U, 3U},
        {2U, 1U, 4U},
        {0U, 2U, 5U}
    };
    mesh.validate();

    const auto coarse_quality = nova::fem::evaluate_mesh_quality_2d(mesh);
    const auto coarse_edges = edge_counts(mesh);
    const auto check = [&](const std::span<const std::size_t> marked,
                           const std::size_t expected_midpoints,
                           const std::size_t expected_children,
                           const std::size_t expected_central_children) {
        const auto refined = nova::fem::refine_triangles_adaptively_2d(
            mesh, marked
        );
        NOVA_REQUIRE(refined.midpoint_parent_edges.size() == expected_midpoints);
        NOVA_REQUIRE(refined.mesh.elements.size() == expected_children);

        std::vector<std::size_t> children_per_parent(mesh.elements.size(), 0U);
        for (const std::size_t parent : refined.parent_elements) {
            NOVA_REQUIRE(parent < mesh.elements.size());
            ++children_per_parent[parent];
        }
        NOVA_REQUIRE(children_per_parent[0] == expected_central_children);

        const auto fine_edges = edge_counts(refined.mesh);
        for (const auto& [edge, incidence] : fine_edges) {
            static_cast<void>(edge);
            NOVA_REQUIRE(incidence == 1U || incidence == 2U);
        }
        for (std::size_t i = 0; i < refined.midpoint_parent_edges.size(); ++i) {
            const BoundaryEdge2D& parent_edge = refined.midpoint_parent_edges[i];
            const std::size_t midpoint = refined.coarse_node_count + i;
            const std::size_t parent_incidence =
                coarse_edges.at(key(parent_edge[0], parent_edge[1]));
            NOVA_REQUIRE(!fine_edges.contains(key(parent_edge[0], parent_edge[1])));
            NOVA_REQUIRE(
                fine_edges.at(key(parent_edge[0], midpoint)) == parent_incidence
            );
            NOVA_REQUIRE(
                fine_edges.at(key(midpoint, parent_edge[1])) == parent_incidence
            );
        }

        const auto fine_quality = nova::fem::evaluate_mesh_quality_2d(refined.mesh);
        NOVA_REQUIRE_NEAR(
            fine_quality.total_area, coarse_quality.total_area, 1.0e-12
        );
        for (std::size_t child = 0; child < refined.parent_elements.size(); ++child) {
            const double parent_area =
                coarse_quality.elements[refined.parent_elements[child]].signed_area;
            NOVA_REQUIRE(
                fine_quality.elements[child].signed_area * parent_area > 0.0
            );
        }
    };

    // Two red-refined neighbors give the central cell two split edges.
    const std::array<std::size_t, 2> two_marks{{1U, 2U}};
    check(two_marks, 6U, 12U, 3U);

    // Three red-refined neighbors give the central cell all three split edges.
    const std::array<std::size_t, 3> three_marks{{1U, 2U, 3U}};
    check(three_marks, 9U, 16U, 4U);
}

NOVA_TEST(adaptive_refinement_preserves_mixed_element_orientation) {
    auto mesh = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    ).mesh;
    std::swap(mesh.elements[0][1], mesh.elements[0][2]);
    const auto coarse = nova::fem::evaluate_mesh_quality_2d(mesh);
    const std::array<std::size_t, 1> marked{{0U}};
    const auto refined = nova::fem::refine_triangles_adaptively_2d(
        mesh, marked
    );
    const auto fine = nova::fem::evaluate_mesh_quality_2d(refined.mesh);

    for (std::size_t child = 0; child < refined.parent_elements.size(); ++child) {
        const double parent_area =
            coarse.elements[refined.parent_elements[child]].signed_area;
        NOVA_REQUIRE(
            fine.elements[child].signed_area * parent_area > 0.0
        );
    }
    NOVA_REQUIRE_NEAR(fine.total_area, coarse.total_area, 1.0e-13);
}

NOVA_TEST(adaptive_refinement_no_marks_returns_identity_transfer) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 2U, 2U
    );
    const auto refined = nova::fem::refine_triangles_adaptively_2d(
        rectangle.mesh, {}, rectangle.boundary_edges
    );
    NOVA_REQUIRE(refined.mesh.nodes.size() == rectangle.mesh.nodes.size());
    NOVA_REQUIRE(refined.mesh.elements.size() == rectangle.mesh.elements.size());
    NOVA_REQUIRE(refined.midpoint_parent_edges.empty());
    NOVA_REQUIRE(refined.parent_elements.size() == rectangle.mesh.elements.size());
    NOVA_REQUIRE(refined.boundary_edges == rectangle.boundary_edges);
}

NOVA_TEST(adaptive_refinement_transfers_quadratic_p2_solution_exactly) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        -1.0, 1.0, -0.5, 1.5, 2U, 2U
    );
    const auto field = [](const double x, const double y) {
        return 0.75 + x - 2.0 * y + x * x + 0.5 * x * y + 2.0 * y * y;
    };
    const H1FunctionSpace2D coarse_space(
        rectangle.mesh, H1TriangleOrder::P2
    );
    const auto coarse_values = coarse_space.interpolate(field);
    const std::array<std::size_t, 2> marked{{1U, 4U}};
    const auto refined = nova::fem::refine_triangles_adaptively_2d(
        rectangle.mesh, marked, rectangle.boundary_edges
    );
    const auto transferred = nova::fem::prolongate_h1_refined_2d(
        rectangle.mesh, refined, H1TriangleOrder::P2, coarse_values
    );
    const H1FunctionSpace2D fine_space(
        refined.mesh, H1TriangleOrder::P2
    );
    const auto exact = fine_space.interpolate(field);
    NOVA_REQUIRE(transferred.size() == exact.size());
    for (std::size_t i = 0; i < exact.size(); ++i) {
        NOVA_REQUIRE_NEAR(transferred[i], exact[i], 3.0e-12);
    }
}

NOVA_TEST(adaptive_thermal_remeshing_preserves_all_dirichlet_edges) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    nova::fem::TransientHeat2DProblem problem;
    problem.mesh = rectangle.mesh;
    problem.volumetric_capacity = [](double, double) { return 1.0; };
    problem.conductivity = [](double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.initial_temperature = [](double, double) { return 0.0; };
    const H1FunctionSpace2D coarse_space(
        problem.mesh, H1TriangleOrder::P2
    );
    problem.initial_dof_values = coarse_space.interpolate(
        [](double, double) { return 0.0; }
    );
    for (const BoundaryEdge2D& edge : rectangle.boundary_edges) {
        problem.dirichlet_edges.push_back({
            edge,
            [](double, double, double) { return 0.0; }
        });
    }

    const std::array<std::size_t, 1> marked{{0U}};
    const auto refined =
        nova::fem::refine_transient_heat_problem_adaptively_2d(
            problem, H1TriangleOrder::P2, marked
        );
    NOVA_REQUIRE(refined.problem.mesh.elements.size() == 6U);
    NOVA_REQUIRE(
        refined.problem.dirichlet_edges.size() ==
        refined.refinement.boundary_edges.size()
    );

    nova::fem::TransientHeat2DOptions options;
    options.order = H1TriangleOrder::P2;
    options.steps = 1U;
    options.time_step = 0.1;
    const auto solution = nova::fem::solve_transient_heat_2d(
        refined.problem, options
    );
    NOVA_REQUIRE(solution.temperature.size() ==
        H1FunctionSpace2D(refined.problem.mesh, H1TriangleOrder::P2).dimension());
}
