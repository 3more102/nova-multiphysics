#pragma once

#include "nova/fem/simplex_quadrature.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace nova::fem::detail {

class ExteriorFaceAssignments3D {
public:
    explicit ExteriorFaceAssignments3D(const TetrahedralMesh3D& mesh)
        : node_count_(mesh.nodes.size()) {
        for (auto face : extract_boundary_faces_3d(mesh)) {
            std::sort(face.begin(), face.end());
            available_.insert(face);
        }
    }

    void claim(BoundaryFace3D face) {
        std::sort(face.begin(), face.end());
        if (face.back() >= node_count_) {
            throw std::out_of_range("Boundary face references a missing vertex.");
        }
        if (!available_.erase(face)) {
            throw std::invalid_argument("Boundary condition requires a distinct, unused exterior face.");
        }
    }

private:
    std::size_t node_count_;
    std::set<BoundaryFace3D> available_;
};

template<std::size_t N>
[[nodiscard]] Point3D simplex_point(
    const TetrahedralMesh3D& mesh, const std::array<std::size_t, N>& nodes,
    const std::array<double, N>& shape
) {
    const auto& origin = mesh.nodes[nodes[0]];
    Point3D offset{};
    for (std::size_t i = 1; i < N; ++i) {
        const auto& p = mesh.nodes[nodes[i]];
        offset.x += (p.x - origin.x) * shape[i];
        offset.y += (p.y - origin.y) * shape[i];
        offset.z += (p.z - origin.z) * shape[i];
    }
    return {origin.x + offset.x, origin.y + offset.y, origin.z + offset.z};
}

[[nodiscard]] inline double face_jacobian(const TetrahedralMesh3D& mesh,
                                           const BoundaryFace3D& face) {
    const auto& a = mesh.nodes[face[0]];
    const auto& b = mesh.nodes[face[1]];
    const auto& c = mesh.nodes[face[2]];
    const double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    const double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    const double jacobian = std::hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
    if (!std::isfinite(jacobian) || jacobian <= 0.0) {
        throw std::invalid_argument("3D boundary face area must be finite and positive.");
    }
    return jacobian;
}

[[nodiscard]] inline double gradient_product(const TetrahedronP1PhysicalGradients3D& basis,
                                              std::size_t i, std::size_t j) {
    const auto& a = basis.gradients[i];
    const auto& b = basis.gradients[j];
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

} // namespace nova::fem::detail
