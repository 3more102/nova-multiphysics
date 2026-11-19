#pragma once

#include "nova/fem/conductivity_3d.hpp"
#include "nova/fem/simplex_quadrature.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <set>
#include <stdexcept>

namespace nova::fem { struct ScalarDiffusion3DProblem; }

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

using P1Matrix3D = std::array<std::array<double, 4>, 4>;

struct ConductivityField3D {
    const SpatialScalarField3D& scalar;
    const SymmetricTensorField3D& tensor;
    std::optional<double> constant{};
    std::optional<SymmetricTensor3D> constant_tensor{};

    void validate() const;
    [[nodiscard]] double scalar_at(const Point3D& point) const;
    [[nodiscard]] Point3D diffusive_flux(const Point3D& point, const Point3D& gradient) const;
};

// Validate the global choice once. The quadrature kernel sees only the material
// of the current cell, including at discontinuities between volume regions.
class DiffusionMaterials3D {
public:
    explicit DiffusionMaterials3D(const ScalarDiffusion3DProblem& problem);
    [[nodiscard]] ConductivityField3D on_cell(std::size_t element) const;

private:
    const ScalarDiffusion3DProblem& problem_;
};

using P1Gradients3D = std::array<std::array<double, 3>, 4>;

void add_conductivity_sample(P1Matrix3D& stiffness, const ConductivityField3D& material,
                             const Point3D& point, const P1Gradients3D& gradients,
                             double sqrt_weight);

// Both physics adapters visit the same physical quadrature points for the
// stiffness and their additional terms (source or consistent mass).
template<class VolumeTerm>
[[nodiscard]] P1Matrix3D integrate_p1_volume_3d(
    const TetrahedralMesh3D& mesh, const Tetrahedron3D& cell,
    const TetrahedronP1PhysicalGradients3D& basis, const ConductivityField3D& material,
    VolumeTerm&& volume_term
) {
    P1Matrix3D stiffness{};
    const double jacobian = std::abs(basis.signed_jacobian_determinant);
    const double sqrt_jacobian = std::sqrt(jacobian);
    for (const auto& q : tetrahedron_degree_two_quadrature) {
        const auto point = simplex_point(mesh, cell, q.shape);
        add_conductivity_sample(stiffness, material, point, basis.gradients,
                                sqrt_jacobian * std::sqrt(q.weight));
        volume_term(point, q.shape, jacobian * q.weight);
    }
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            stiffness[j][i] = stiffness[i][j];
        }
    }
    return stiffness;
}

} // namespace nova::fem::detail
