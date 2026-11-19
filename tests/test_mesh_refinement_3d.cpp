#include "test_framework.hpp"
#include "nova/fem/mesh_refinement_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

[[nodiscard]] double volume(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::Tetrahedron3D& t
) {
    const auto& p = mesh.nodes[t[0]];
    const auto& q = mesh.nodes[t[1]];
    const auto& r = mesh.nodes[t[2]];
    const auto& s = mesh.nodes[t[3]];
    const double ax = q.x - p.x, ay = q.y - p.y, az = q.z - p.z;
    const double bx = r.x - p.x, by = r.y - p.y, bz = r.z - p.z;
    const double cx = s.x - p.x, cy = s.y - p.y, cz = s.z - p.z;
    return std::abs(ax * (by * cz - bz * cy) -
                    ay * (bx * cz - bz * cx) +
                    az * (bx * cy - by * cx)) / 6.0;
}

[[nodiscard]] std::array<std::size_t, 3> sorted_face(nova::fem::BoundaryFace3D f) {
    std::sort(f.begin(), f.end());
    return f;
}

[[nodiscard]] double outward_dot(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::BoundaryFace3D& face,
    const nova::fem::Point3D& center
) {
    const auto& a = mesh.nodes[face[0]];
    const auto& b = mesh.nodes[face[1]];
    const auto& c = mesh.nodes[face[2]];
    const double x1 = b.x - a.x, y1 = b.y - a.y, z1 = b.z - a.z;
    const double x2 = c.x - a.x, y2 = c.y - a.y, z2 = c.z - a.z;
    const double nx = y1 * z2 - z1 * y2;
    const double ny = z1 * x2 - x1 * z2;
    const double nz = x1 * y2 - y1 * x2;
    return nx * (a.x - center.x) + ny * (a.y - center.y) +
        nz * (a.z - center.z);
}

} // namespace

NOVA_TEST(uniform_tetra_refinement_preserves_volume_boundary_and_parent_maps) {
    const nova::fem::TetrahedralMesh3D mesh{
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
        {{0U, 1U, 2U, 3U}}
    };
    const auto boundary = nova::fem::extract_boundary_faces_3d(mesh);
    const auto result = nova::fem::refine_tetrahedra_uniformly_3d(mesh, boundary);
    NOVA_REQUIRE(result.coarse_node_count == 4U);
    NOVA_REQUIRE(result.coarse_element_count == 1U);
    NOVA_REQUIRE(result.mesh.nodes.size() == 10U);
    NOVA_REQUIRE(result.midpoint_parent_edges.size() == 6U);
    NOVA_REQUIRE(result.mesh.elements.size() == 8U);
    NOVA_REQUIRE(result.parent_elements.size() == 8U);
    NOVA_REQUIRE(result.boundary_faces.size() == 16U);
    NOVA_REQUIRE(result.parent_boundary_faces.size() == 16U);
    result.mesh.validate();

    double refined_volume = 0.0;
    for (std::size_t i = 0; i < result.mesh.elements.size(); ++i) {
        const double child_volume = volume(result.mesh, result.mesh.elements[i]);
        NOVA_REQUIRE(child_volume > 0.0);
        NOVA_REQUIRE_NEAR(child_volume, 1.0 / 48.0, 1.0e-14);
        refined_volume += child_volume;
        NOVA_REQUIRE(result.parent_elements[i] == 0U);
    }
    NOVA_REQUIRE_NEAR(refined_volume, volume(mesh, mesh.elements[0]), 1.0e-14);

    const auto extracted = nova::fem::extract_boundary_faces_3d(result.mesh);
    std::vector<std::array<std::size_t, 3>> refined_keys, extracted_keys;
    for (const auto& face : result.boundary_faces) {
        NOVA_REQUIRE(outward_dot(result.mesh, face, {0.25, 0.25, 0.25}) > 0.0);
        refined_keys.push_back(sorted_face(face));
    }
    for (const auto& face : extracted) {
        extracted_keys.push_back(sorted_face(face));
    }
    std::sort(refined_keys.begin(), refined_keys.end());
    std::sort(extracted_keys.begin(), extracted_keys.end());
    NOVA_REQUIRE(refined_keys == extracted_keys);
    for (std::size_t i = 0; i < result.parent_boundary_faces.size(); ++i) {
        NOVA_REQUIRE(result.parent_boundary_faces[i] == i / 4U);
    }
}

NOVA_TEST(uniform_tetra_refinement_conforms_across_cells_and_two_levels) {
    const nova::fem::TetrahedralMesh3D mesh{
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}},
        {{0U, 1U, 2U, 3U}, {0U, 2U, 1U, 4U}}
    };
    const auto first = nova::fem::refine_tetrahedra_uniformly_3d(mesh);
    NOVA_REQUIRE(first.mesh.elements.size() == 16U);
    NOVA_REQUIRE(first.mesh.nodes.size() == 14U);
    double first_volume = 0.0;
    for (std::size_t i = 0; i < first.mesh.elements.size(); ++i) {
        first_volume += volume(first.mesh, first.mesh.elements[i]);
        NOVA_REQUIRE(first.parent_elements[i] == i / 8U);
    }
    NOVA_REQUIRE_NEAR(first_volume, 1.0 / 3.0, 1.0e-13);
    const auto second = nova::fem::refine_tetrahedra_uniformly_3d(first.mesh);
    NOVA_REQUIRE(second.mesh.elements.size() == 128U);
    double second_volume = 0.0;
    for (const auto& tet : second.mesh.elements) {
        second_volume += volume(second.mesh, tet);
    }
    NOVA_REQUIRE_NEAR(second_volume, 1.0 / 3.0, 1.0e-13);
    second.mesh.validate();
}

NOVA_TEST(uniform_tetra_refinement_transfers_affine_p1_and_cell_data) {
    const nova::fem::TetrahedralMesh3D mesh{
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}},
        {{0U, 1U, 2U, 3U}, {0U, 2U, 1U, 4U}}
    };
    const auto refined = nova::fem::refine_tetrahedra_uniformly_3d(mesh);
    std::vector<double> coarse;
    for (const auto& p : mesh.nodes) {
        coarse.push_back(2.0 + 3.0 * p.x - 4.0 * p.y + 5.0 * p.z);
    }
    const auto values = nova::fem::prolongate_p1_uniform_3d(refined, coarse);
    NOVA_REQUIRE(values.size() == refined.mesh.nodes.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        const auto& p = refined.mesh.nodes[i];
        NOVA_REQUIRE_NEAR(values[i], 2.0 + 3.0 * p.x - 4.0 * p.y + 5.0 * p.z, 1e-13);
    }
    const auto cells = nova::fem::prolongate_cellwise_constant_uniform_3d(
        refined, std::vector<double>{2.5, 7.0}
    );
    NOVA_REQUIRE(cells.size() == 16U);
    for (std::size_t i = 0; i < cells.size(); ++i) {
        NOVA_REQUIRE_NEAR(cells[i], i < 8U ? 2.5 : 7.0, 0.0);
    }

    auto shuffled = mesh;
    std::reverse(shuffled.elements.begin(), shuffled.elements.end());
    std::swap(shuffled.elements[0][2], shuffled.elements[0][3]);
    const auto reordered = nova::fem::refine_tetrahedra_uniformly_3d(shuffled);
    NOVA_REQUIRE(reordered.midpoint_parent_edges == refined.midpoint_parent_edges);
    NOVA_REQUIRE(reordered.mesh.nodes.size() == refined.mesh.nodes.size());
    for (std::size_t i = 0; i < refined.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(reordered.mesh.nodes[i].x, refined.mesh.nodes[i].x, 0.0);
        NOVA_REQUIRE_NEAR(reordered.mesh.nodes[i].y, refined.mesh.nodes[i].y, 0.0);
        NOVA_REQUIRE_NEAR(reordered.mesh.nodes[i].z, refined.mesh.nodes[i].z, 0.0);
    }
}

NOVA_TEST(uniform_tetra_refinement_rejects_bad_boundaries_and_fields) {
    const nova::fem::TetrahedralMesh3D mesh{
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
        {{0U, 1U, 2U, 3U}}
    };
    const auto faces = nova::fem::extract_boundary_faces_3d(mesh);
    auto reversed = faces.front();
    std::swap(reversed[0], reversed[1]);
    bool rejected = false;
    try {
        static_cast<void>(nova::fem::refine_tetrahedra_uniformly_3d(mesh, {reversed}));
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    rejected = false;
    try {
        static_cast<void>(nova::fem::refine_tetrahedra_uniformly_3d(mesh, {faces[0], faces[0]}));
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    auto bad_face = faces[0];
    bad_face[1] = mesh.nodes.size();
    rejected = false;
    try {
        static_cast<void>(nova::fem::refine_tetrahedra_uniformly_3d(mesh, {bad_face}));
    } catch (const std::out_of_range&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    const auto refined = nova::fem::refine_tetrahedra_uniformly_3d(mesh);
    rejected = false;
    try {
        static_cast<void>(nova::fem::prolongate_p1_uniform_3d(
            refined, std::vector<double>{1.0}
        ));
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    auto coarse = std::vector<double>(4U, 2.0);
    coarse[1] = std::numeric_limits<double>::quiet_NaN();
    rejected = false;
    try {
        static_cast<void>(nova::fem::prolongate_p1_uniform_3d(refined, coarse));
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);

    rejected = false;
    try {
        static_cast<void>(nova::fem::prolongate_cellwise_constant_uniform_3d(
            refined, std::vector<double>{}
        ));
    } catch (const std::invalid_argument&) { rejected = true; }
    NOVA_REQUIRE(rejected);
}
