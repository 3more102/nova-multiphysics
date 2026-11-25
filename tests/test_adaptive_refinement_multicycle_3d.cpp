#include "test_framework.hpp"

#include "nova/fem/mesh_refinement_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

namespace {

using nova::fem::BoundaryFace3D;
using nova::fem::RefinedTetrahedralMesh3D;
using FaceKey = std::array<std::size_t, 3>;

[[nodiscard]] FaceKey canonical(FaceKey face) {
    std::sort(face.begin(), face.end());
    return face;
}

[[nodiscard]] bool same_outward_winding(
    const BoundaryFace3D& actual, const BoundaryFace3D& expected
) {
    return actual == expected ||
        (actual[0] == expected[1] && actual[1] == expected[2] &&
         actual[2] == expected[0]) ||
        (actual[0] == expected[2] && actual[1] == expected[0] &&
         actual[2] == expected[1]);
}

// Detect missing/duplicate exterior faces as well as incorrect orientation.
void verify_boundary_partition(
    const RefinedTetrahedralMesh3D& stage,
    const std::size_t parent_face_count
) {
    const auto actual_exterior =
        nova::fem::extract_boundary_faces_3d(stage.mesh);
    NOVA_REQUIRE(stage.boundary_faces.size() == actual_exterior.size());
    NOVA_REQUIRE(stage.parent_boundary_faces.size() == stage.boundary_faces.size());
    std::map<FaceKey, BoundaryFace3D> exterior;
    for (const auto& face : actual_exterior) {
        NOVA_REQUIRE(exterior.emplace(canonical(face), face).second);
    }
    std::map<FaceKey, bool> seen;
    for (std::size_t i = 0; i < stage.boundary_faces.size(); ++i) {
        const auto& face = stage.boundary_faces[i];
        const auto found = exterior.find(canonical(face));
        NOVA_REQUIRE(found != exterior.end());
        NOVA_REQUIRE(seen.emplace(canonical(face), true).second);
        NOVA_REQUIRE(same_outward_winding(face, found->second));
        NOVA_REQUIRE(stage.parent_boundary_faces[i] < parent_face_count);
    }
}

} // namespace

NOVA_TEST(tetra_two_local_refinements_preserve_material_volume_and_affine_field) {
    // The second refinement is performed on a mesh that is already locally
    // subdivided; single-step checks alone cannot detect mapping drift.
    const auto original = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const auto first = nova::fem::refine_tetrahedra_adaptively_3d(
        original.mesh, std::vector<std::size_t>{0U}, original.boundary_faces);
    first.mesh.validate();
    NOVA_REQUIRE(first.mesh.elements.size() > original.mesh.elements.size());

    const auto second = nova::fem::refine_tetrahedra_adaptively_3d(
        first.mesh, std::vector<std::size_t>{first.mesh.elements.size() / 2U},
        first.boundary_faces);
    second.mesh.validate();
    NOVA_REQUIRE(second.mesh.elements.size() > first.mesh.elements.size());

    verify_boundary_partition(first, original.boundary_faces.size());
    verify_boundary_partition(second, first.boundary_faces.size());

    // Trace every final exterior triangle through both parent maps to a
    // legitimate original boundary face.
    for (const auto parent1 : second.parent_boundary_faces) {
        NOVA_REQUIRE(parent1 < first.parent_boundary_faces.size());
        const auto parent0 = first.parent_boundary_faces[parent1];
        NOVA_REQUIRE(parent0 < original.boundary_faces.size());
    }

    const auto affine = [](double x, double y, double z) {
        return 1.75 + 2.0 * x - 3.0 * y + 0.5 * z;
    };
    std::vector<double> coarse_nodal;
    for (const auto& p : original.mesh.nodes) {
        coarse_nodal.push_back(affine(p.x, p.y, p.z));
    }
    const auto first_nodal =
        nova::fem::prolongate_p1_refined_3d(first, coarse_nodal);
    const auto second_nodal =
        nova::fem::prolongate_p1_refined_3d(second, first_nodal);
    NOVA_REQUIRE(second_nodal.size() == second.mesh.nodes.size());
    for (std::size_t i = 0; i < second_nodal.size(); ++i) {
        const auto& p = second.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(
            second_nodal[i], affine(p.x, p.y, p.z), 1.0e-12);
    }

    std::vector<double> material;
    material.reserve(original.mesh.elements.size());
    for (std::size_t i = 0; i < original.mesh.elements.size(); ++i) {
        material.push_back(1.0 + 0.375 * static_cast<double>(i));
    }
    const auto first_material =
        nova::fem::prolongate_cellwise_constant_refined_3d(first, material);
    const auto second_material =
        nova::fem::prolongate_cellwise_constant_refined_3d(second, first_material);
    NOVA_REQUIRE(second_material.size() == second.mesh.elements.size());

    const auto coarse_quality = nova::fem::evaluate_mesh_quality_3d(original.mesh);
    const auto fine_quality = nova::fem::evaluate_mesh_quality_3d(second.mesh);
    std::vector<long double> volumes(original.mesh.elements.size(), 0.0L);
    std::vector<long double> weighted(original.mesh.elements.size(), 0.0L);
    for (std::size_t i = 0; i < second.mesh.elements.size(); ++i) {
        const auto immediate_parent = second.parent_elements[i];
        NOVA_REQUIRE(immediate_parent < first.parent_elements.size());
        const auto root = first.parent_elements[immediate_parent];
        NOVA_REQUIRE(root < material.size());
        NOVA_REQUIRE(second_material[i] == material[root]);
        const auto child_volume =
            static_cast<long double>(fine_quality.elements[i].volume);
        volumes[root] += child_volume;
        weighted[root] += child_volume * second_material[i];
    }
    for (std::size_t root = 0; root < material.size(); ++root) {
        const auto expected_volume =
            static_cast<long double>(coarse_quality.elements[root].volume);
        NOVA_REQUIRE_NEAR(
            static_cast<double>(volumes[root]),
            static_cast<double>(expected_volume), 1.0e-12);
        NOVA_REQUIRE_NEAR(
            static_cast<double>(weighted[root]),
            static_cast<double>(expected_volume * material[root]), 1.0e-12);
    }

    // Both steps are deterministic for identical meshes and marking sets.
    const auto second_again = nova::fem::refine_tetrahedra_adaptively_3d(
        first.mesh, std::vector<std::size_t>{first.mesh.elements.size() / 2U},
        first.boundary_faces);
    NOVA_REQUIRE(second_again.mesh.elements == second.mesh.elements);
    NOVA_REQUIRE(second_again.parent_elements == second.parent_elements);
    NOVA_REQUIRE(second_again.midpoint_parent_edges == second.midpoint_parent_edges);
    NOVA_REQUIRE(second_again.boundary_faces == second.boundary_faces);
    NOVA_REQUIRE(second_again.parent_boundary_faces ==
                 second.parent_boundary_faces);
}
