#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace nova::fem {

struct Point3D {
    double x{};
    double y{};
    double z{};
};

using Tetrahedron3D = std::array<std::size_t, 4>;
using BoundaryFace3D = std::array<std::size_t, 3>;

struct TetrahedralMesh3D {
    std::vector<Point3D> nodes;
    std::vector<Tetrahedron3D> elements;

    // Checks finite coordinates, valid unique connectivity, non-degenerate
    // tetrahedra, and manifold face incidence.
    void validate() const;
};

struct StructuredBoxMesh3D {
    TetrahedralMesh3D mesh;
    // Exterior faces are oriented outward from the mesh domain.
    std::vector<BoundaryFace3D> boundary_faces;
};

struct TetrahedronQuality3D {
    double signed_volume{};
    double volume{};
    double mean_ratio{};
    double minimum_dihedral_angle{}; // radians
    double maximum_edge_length{};
};

struct MeshQuality3D {
    std::vector<TetrahedronQuality3D> elements;
    double total_volume{};
    double minimum_mean_ratio{};
    double minimum_dihedral_angle{}; // radians
};

// Creates a structured box split into six conforming P1 tetrahedra per cell.
[[nodiscard]] StructuredBoxMesh3D make_box_mesh_3d(
    double x_min, double x_max,
    double y_min, double y_max,
    double z_min, double z_max,
    std::size_t elements_x,
    std::size_t elements_y,
    std::size_t elements_z
);

// Returns each exterior triangular face exactly once, with outward orientation.
[[nodiscard]] std::vector<BoundaryFace3D> extract_boundary_faces_3d(
    const TetrahedralMesh3D& mesh
);

[[nodiscard]] MeshQuality3D evaluate_mesh_quality_3d(
    const TetrahedralMesh3D& mesh
);

} // namespace nova::fem
