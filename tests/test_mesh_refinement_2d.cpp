#include "test_framework.hpp"

#include "nova/fem/mesh_refinement_2d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using nova::fem::BoundaryEdge2D;
using nova::fem::TriangleMesh2D;

[[nodiscard]] std::pair<std::size_t, std::size_t> key(
    std::size_t first, std::size_t second
) {
    return std::minmax(first, second);
}

[[nodiscard]] std::map<std::pair<std::size_t, std::size_t>, std::size_t>
edge_counts(const TriangleMesh2D& mesh) {
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> counts;
    for (const auto& cell : mesh.elements) {
        ++counts[key(cell[0], cell[1])];
        ++counts[key(cell[1], cell[2])];
        ++counts[key(cell[2], cell[0])];
    }
    return counts;
}
} // namespace

NOVA_TEST(mesh_quality_equilateral_and_right_triangle) {
    constexpr double pi = 3.14159265358979323846;
    const TriangleMesh2D equilateral{{
        {0.0, 0.0}, {1.0, 0.0}, {0.5, std::sqrt(3.0) / 2.0}
    }, {{0U, 1U, 2U}}};
    const auto equi = nova::fem::evaluate_mesh_quality_2d(equilateral);
    NOVA_REQUIRE(equi.elements.size() == 1U);
    NOVA_REQUIRE_NEAR(equi.minimum_mean_ratio_quality, 1.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(equi.minimum_angle_radians, pi / 3.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(equi.total_area, std::sqrt(3.0) / 4.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(equi.maximum_edge_length, 1.0, 1.0e-14);

    const TriangleMesh2D right{{
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}
    }, {{0U, 1U, 2U}}};
    const auto q = nova::fem::evaluate_mesh_quality_2d(right);
    NOVA_REQUIRE_NEAR(q.minimum_mean_ratio_quality, std::sqrt(3.0) / 2.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(q.minimum_angle_radians, pi / 4.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(q.total_area, 0.5, 1.0e-14);
}

NOVA_TEST(mesh_uniform_refinement_preserves_conforming_topology_and_area) {
    const auto original = nova::fem::make_rectangle_mesh_2d(
        -1.0, 2.0, 0.0, 1.0, 2U, 3U
    );
    const auto coarse_edges = edge_counts(original.mesh);
    const auto refined = nova::fem::refine_triangles_uniformly_2d(
        original.mesh, original.boundary_edges
    );
    NOVA_REQUIRE(refined.coarse_node_count == original.mesh.nodes.size());
    NOVA_REQUIRE(refined.mesh.nodes.size() ==
        original.mesh.nodes.size() + coarse_edges.size());
    NOVA_REQUIRE(refined.midpoint_parent_edges.size() == coarse_edges.size());
    NOVA_REQUIRE(refined.mesh.elements.size() ==
        4U * original.mesh.elements.size());
    NOVA_REQUIRE(refined.parent_elements.size() ==
        refined.mesh.elements.size());
    NOVA_REQUIRE(refined.boundary_edges.size() ==
        2U * original.boundary_edges.size());

    const auto old_quality = nova::fem::evaluate_mesh_quality_2d(original.mesh);
    const auto new_quality = nova::fem::evaluate_mesh_quality_2d(refined.mesh);
    NOVA_REQUIRE_NEAR(old_quality.total_area, new_quality.total_area, 1.0e-12);
    NOVA_REQUIRE_NEAR(
        old_quality.minimum_mean_ratio_quality,
        new_quality.minimum_mean_ratio_quality, 1.0e-12
    );
    NOVA_REQUIRE_NEAR(
        old_quality.minimum_angle_radians,
        new_quality.minimum_angle_radians, 1.0e-12
    );
    for (std::size_t i = 0; i < refined.parent_elements.size(); ++i) {
        NOVA_REQUIRE(refined.parent_elements[i] == i / 4U);
    }

    const auto refined_edges = edge_counts(refined.mesh);
    for (std::size_t i = 0; i < original.boundary_edges.size(); ++i) {
        const BoundaryEdge2D coarse = original.boundary_edges[i];
        const auto left = refined.boundary_edges[2U * i];
        const auto right = refined.boundary_edges[2U * i + 1U];
        NOVA_REQUIRE(left[0] == coarse[0]);
        NOVA_REQUIRE(left[1] == right[0]);
        NOVA_REQUIRE(right[1] == coarse[1]);
        NOVA_REQUIRE(refined_edges.at(key(left[0], left[1])) == 1U);
        NOVA_REQUIRE(refined_edges.at(key(right[0], right[1])) == 1U);
    }
    for (const auto& [edge, incidence] : refined_edges) {
        static_cast<void>(edge);
        NOVA_REQUIRE(incidence == 1U || incidence == 2U);
    }
}

NOVA_TEST(mesh_uniform_refinement_prolongates_affine_p1_exactly) {
    const auto original = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 2.0, 2U, 2U
    );
    const auto refined = nova::fem::refine_triangles_uniformly_2d(
        original.mesh
    );
    std::vector<double> coarse;
    for (const auto& node : original.mesh.nodes) {
        coarse.push_back(1.5 + 2.0 * node.x - 3.0 * node.y);
    }
    const auto fine = nova::fem::prolongate_p1_uniform_2d(refined, coarse);
    NOVA_REQUIRE(fine.size() == refined.mesh.nodes.size());
    for (std::size_t i = 0; i < fine.size(); ++i) {
        const auto& node = refined.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(
            fine[i], 1.5 + 2.0 * node.x - 3.0 * node.y, 1.0e-13
        );
    }
    NOVA_REQUIRE(refined.boundary_edges.empty());
}

NOVA_TEST(mesh_uniform_refinement_retains_clockwise_orientation) {
    auto mesh = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    ).mesh;
    std::swap(mesh.elements[0][1], mesh.elements[0][2]);
    const auto coarse = nova::fem::evaluate_mesh_quality_2d(mesh);
    NOVA_REQUIRE(coarse.elements[0].signed_area < 0.0);
    NOVA_REQUIRE(coarse.elements[1].signed_area > 0.0);
    const auto refined = nova::fem::refine_triangles_uniformly_2d(mesh);
    const auto fine = nova::fem::evaluate_mesh_quality_2d(refined.mesh);
    for (std::size_t i = 0; i < 4U; ++i) {
        NOVA_REQUIRE(fine.elements[i].signed_area < 0.0);
    }
    for (std::size_t i = 4U; i < 8U; ++i) {
        NOVA_REQUIRE(fine.elements[i].signed_area > 0.0);
    }
    NOVA_REQUIRE_NEAR(fine.total_area, coarse.total_area, 1.0e-14);
}

NOVA_TEST(mesh_uniform_refinement_rejects_interior_or_duplicate_boundaries) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    const auto counts = edge_counts(rectangle.mesh);
    std::pair<std::size_t, std::size_t> interior{};
    for (const auto& [edge, incidence] : counts) {
        if (incidence == 2U) {
            interior = edge;
        }
    }
    bool rejected_interior = false;
    try {
        static_cast<void>(nova::fem::refine_triangles_uniformly_2d(
            rectangle.mesh, {{interior.first, interior.second}}
        ));
    } catch (const std::invalid_argument&) {
        rejected_interior = true;
    }
    NOVA_REQUIRE(rejected_interior);

    const auto edge = rectangle.boundary_edges.front();
    bool rejected_duplicate = false;
    try {
        static_cast<void>(nova::fem::refine_triangles_uniformly_2d(
            rectangle.mesh, {edge, {edge[1], edge[0]}}
        ));
    } catch (const std::invalid_argument&) {
        rejected_duplicate = true;
    }
    NOVA_REQUIRE(rejected_duplicate);
}

NOVA_TEST(mesh_uniform_refinement_rejects_bad_p1_state) {
    const auto mesh = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    const auto refined = nova::fem::refine_triangles_uniformly_2d(mesh.mesh);
    bool rejected_size = false;
    try {
        static_cast<void>(nova::fem::prolongate_p1_uniform_2d(
            refined, std::vector<double>(1U, 1.0)
        ));
    } catch (const std::invalid_argument&) {
        rejected_size = true;
    }
    NOVA_REQUIRE(rejected_size);
    std::vector<double> values(mesh.mesh.nodes.size(), 1.0);
    values[0] = std::numeric_limits<double>::quiet_NaN();
    bool rejected_nan = false;
    try {
        static_cast<void>(nova::fem::prolongate_p1_uniform_2d(refined, values));
    } catch (const std::invalid_argument&) {
        rejected_nan = true;
    }
    NOVA_REQUIRE(rejected_nan);
}
