#include "test_framework.hpp"
#include "nova/fem/mesh_refinement_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

namespace {
using namespace nova::fem;

template<class F>
bool invalid(F&& fn) {
    try { fn(); } catch (const std::invalid_argument&) { return true; }
    return false;
}

[[nodiscard]] BoundaryFace3D sorted(BoundaryFace3D f) {
    std::sort(f.begin(), f.end());
    return f;
}

[[nodiscard]] double affine(const Point3D& p) {
    return 2.0 + 3.0 * p.x - 4.0 * p.y + 5.0 * p.z;
}

} // namespace

NOVA_TEST(tetrahedron_uniform_refinement_has_eight_equal_volume_children) {
    const TetrahedralMesh3D original{
        {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}},
        {{0,1,2,3}}
    };
    for (bool reverse : {false, true}) {
        auto coarse = original;
        if (reverse) std::swap(coarse.elements[0][1], coarse.elements[0][2]);
        const auto refined = refine_tetrahedra_uniform_3d(coarse);
        NOVA_REQUIRE(refined.original_node_count == 4U);
        NOVA_REQUIRE(refined.mesh.nodes.size() == 10U);
        NOVA_REQUIRE(refined.midpoint_parent_edges.size() == 6U);
        NOVA_REQUIRE(refined.mesh.elements.size() == 8U);
        NOVA_REQUIRE(refined.parent_elements.size() == 8U);
        NOVA_REQUIRE(refined.boundary_faces.size() == 16U);
        NOVA_REQUIRE(refined.parent_boundary_faces.size() == 16U);
        for (std::size_t i = 0; i < 8U; ++i) {
            NOVA_REQUIRE(refined.parent_elements[i] == 0U);
        }
        const auto coarse_q = evaluate_mesh_quality_3d(coarse);
        const auto fine_q = evaluate_mesh_quality_3d(refined.mesh);
        NOVA_REQUIRE_NEAR(fine_q.total_volume, coarse_q.total_volume, 1e-14);
        for (const auto& quality : fine_q.elements) {
            NOVA_REQUIRE_NEAR(quality.volume, 1.0 / 48.0, 1e-14);
            NOVA_REQUIRE(quality.mean_ratio > 0.0);
        }
        for (std::size_t i = 0; i < 16U; ++i) {
            NOVA_REQUIRE(refined.parent_boundary_faces[i] == i / 4U);
        }
        const auto extracted = extract_boundary_faces_3d(refined.mesh);
        std::multiset<BoundaryFace3D> actual, mapped;
        for (const auto f : extracted) actual.insert(sorted(f));
        for (const auto f : refined.boundary_faces) mapped.insert(sorted(f));
        NOVA_REQUIRE(actual == mapped);
    }
}

NOVA_TEST(tetrahedron_uniform_refinement_shares_midpoints_across_cells) {
    // Two tetrahedra occupy opposite sides of the same triangular face.
    const TetrahedralMesh3D coarse{
        {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {0,0,-1}},
        {{0,1,2,3}, {0,2,1,4}}
    };
    const auto refined = refine_tetrahedra_uniform_3d(coarse);
    NOVA_REQUIRE(refined.mesh.nodes.size() == 14U); // 5 vertices, 9 unique edges
    NOVA_REQUIRE(refined.mesh.elements.size() == 16U);
    NOVA_REQUIRE(refined.boundary_faces.size() == 24U);
    const auto before = evaluate_mesh_quality_3d(coarse);
    const auto after = evaluate_mesh_quality_3d(refined.mesh);
    NOVA_REQUIRE_NEAR(after.total_volume, before.total_volume, 1e-13);
    std::multiset<BoundaryFace3D> actual, mapped;
    for (const auto f : extract_boundary_faces_3d(refined.mesh)) actual.insert(sorted(f));
    for (const auto f : refined.boundary_faces) mapped.insert(sorted(f));
    NOVA_REQUIRE(actual == mapped);
    NOVA_REQUIRE(refined.parent_elements[0] == 0U);
    NOVA_REQUIRE(refined.parent_elements[8] == 1U);
}

NOVA_TEST(tetrahedron_p1_transfer_is_exact_for_affine_fields) {
    const auto coarse = make_box_mesh_3d(0, 2, -1, 1, 0, 3, 2, 1, 1).mesh;
    const auto refined = refine_tetrahedra_uniform_3d(coarse);
    std::vector<double> field;
    for (const auto& node : coarse.nodes) field.push_back(affine(node));
    const auto fine = prolongate_tetrahedral_p1_3d(refined, field);
    NOVA_REQUIRE(fine.size() == refined.mesh.nodes.size());
    for (std::size_t i = 0; i < fine.size(); ++i) {
        NOVA_REQUIRE_NEAR(fine[i], affine(refined.mesh.nodes[i]), 1e-13);
    }
}

NOVA_TEST(tetrahedron_p1_transfer_rejects_invalid_metadata_and_values) {
    const auto coarse = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 1, 1, 1).mesh;
    auto refined = refine_tetrahedra_uniform_3d(coarse);
    std::vector<double> input(coarse.nodes.size(), 1.0);
    NOVA_REQUIRE(invalid([&] {
        auto short_field = input; short_field.pop_back();
        (void)prolongate_tetrahedral_p1_3d(refined, short_field);
    }));
    NOVA_REQUIRE(invalid([&] {
        auto bad = input;
        bad[0] = std::numeric_limits<double>::quiet_NaN();
        (void)prolongate_tetrahedral_p1_3d(refined, bad);
    }));
    NOVA_REQUIRE(invalid([&] {
        auto damaged = refined; damaged.midpoint_parent_edges[0][0] = input.size();
        (void)prolongate_tetrahedral_p1_3d(damaged, input);
    }));
    NOVA_REQUIRE(invalid([&] {
        auto damaged = refined; damaged.original_node_count = input.size() + 1U;
        (void)prolongate_tetrahedral_p1_3d(damaged, input);
    }));
}
