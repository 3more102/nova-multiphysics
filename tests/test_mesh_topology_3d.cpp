#include "test_framework.hpp"

#include "nova/fem/mesh3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

[[nodiscard]] double face_normal_dot_outward(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::BoundaryFace3D& face,
    const nova::fem::Point3D& center
) {
    const auto& a = mesh.nodes[face[0]];
    const auto& b = mesh.nodes[face[1]];
    const auto& c = mesh.nodes[face[2]];
    const double ux = b.x - a.x;
    const double uy = b.y - a.y;
    const double uz = b.z - a.z;
    const double vx = c.x - a.x;
    const double vy = c.y - a.y;
    const double vz = c.z - a.z;
    const double nx = uy * vz - uz * vy;
    const double ny = uz * vx - ux * vz;
    const double nz = ux * vy - uy * vx;
    const double cx = (a.x + b.x + c.x) / 3.0;
    const double cy = (a.y + b.y + c.y) / 3.0;
    const double cz = (a.z + b.z + c.z) / 3.0;
    return nx * (cx - center.x) +
        ny * (cy - center.y) +
        nz * (cz - center.z);
}

} // namespace

NOVA_TEST(box_mesh_3d_has_six_tets_per_cell_and_outward_boundary_faces) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U
    );
    NOVA_REQUIRE(box.mesh.nodes.size() == 8U);
    NOVA_REQUIRE(box.mesh.elements.size() == 6U);
    NOVA_REQUIRE(box.boundary_faces.size() == 12U);

    const auto quality = nova::fem::evaluate_mesh_quality_3d(box.mesh);
    NOVA_REQUIRE_NEAR(quality.total_volume, 1.0, 1.0e-13);
    NOVA_REQUIRE(quality.minimum_mean_ratio > 0.0);
    NOVA_REQUIRE(quality.minimum_dihedral_angle > 0.0);

    const nova::fem::Point3D center{0.5, 0.5, 0.5};
    for (const auto& face : box.boundary_faces) {
        NOVA_REQUIRE(face_normal_dot_outward(box.mesh, face, center) > 0.0);
    }
}

NOVA_TEST(box_mesh_3d_is_conforming_across_multiple_cells) {
    const auto box = nova::fem::make_box_mesh_3d(
        -1.0, 2.0, 0.0, 2.0, -0.5, 1.5, 2U, 2U, 2U
    );
    NOVA_REQUIRE(box.mesh.elements.size() == 48U);
    NOVA_REQUIRE_NEAR(
        nova::fem::evaluate_mesh_quality_3d(box.mesh).total_volume,
        12.0, 1.0e-12
    );
    NOVA_REQUIRE(box.boundary_faces.size() == 48U);
}

NOVA_TEST(tetrahedral_mesh_3d_accepts_both_orientations) {
    auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U
    );
    std::swap(box.mesh.elements[0][2], box.mesh.elements[0][3]);
    box.mesh.validate();
    NOVA_REQUIRE_NEAR(
        nova::fem::evaluate_mesh_quality_3d(box.mesh).total_volume,
        1.0, 1.0e-13
    );
}

NOVA_TEST(regular_tetrahedron_has_unit_mean_ratio_and_expected_dihedral) {
    const double root3 = std::sqrt(3.0);
    const double root23 = std::sqrt(2.0 / 3.0);
    const nova::fem::TetrahedralMesh3D mesh{
        {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.5, root3 / 2.0, 0.0},
            {0.5, root3 / 6.0, root23}
        },
        {{0U, 1U, 2U, 3U}}
    };
    const auto quality = nova::fem::evaluate_mesh_quality_3d(mesh);
    NOVA_REQUIRE_NEAR(quality.minimum_mean_ratio, 1.0, 2.0e-13);
    NOVA_REQUIRE_NEAR(
        quality.minimum_dihedral_angle,
        1.2309594173407747,
        2.0e-12
    );
}

NOVA_TEST(tetrahedral_mesh_3d_rejects_degenerate_and_nonmanifold_cells) {
    const nova::fem::TetrahedralMesh3D degenerate{
        {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {1.0, 1.0, 0.0}
        },
        {{0U, 1U, 2U, 3U}}
    };
    bool rejected_degenerate = false;
    try {
        degenerate.validate();
    } catch (const std::invalid_argument&) {
        rejected_degenerate = true;
    }
    NOVA_REQUIRE(rejected_degenerate);

    const nova::fem::TetrahedralMesh3D nonmanifold{
        {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, -1.0},
            {1.0, 1.0, 1.0}
        },
        {
            {0U, 1U, 2U, 3U},
            {0U, 2U, 1U, 4U},
            {0U, 1U, 2U, 5U}
        }
    };
    bool rejected_nonmanifold = false;
    try {
        nonmanifold.validate();
    } catch (const std::invalid_argument&) {
        rejected_nonmanifold = true;
    }
    NOVA_REQUIRE(rejected_nonmanifold);
}
