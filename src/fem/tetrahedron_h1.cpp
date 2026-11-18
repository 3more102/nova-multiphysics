#include "nova/fem/tetrahedron_h1.hpp"

#include <cmath>
#include <stdexcept>

namespace nova::fem {
namespace {

struct Vector3 {
    double x{};
    double y{};
    double z{};
};

[[nodiscard]] Vector3 subtract(const Point3D& a, const Point3D& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

[[nodiscard]] Vector3 cross(const Vector3& a, const Vector3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

[[nodiscard]] double dot(const Vector3& a, const Vector3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

} // namespace

TetrahedronP1BasisEvaluation3D evaluate_tetrahedron_p1_basis_3d(
    const double r,
    const double s,
    const double t
) {
    if (!std::isfinite(r) || !std::isfinite(s) || !std::isfinite(t)) {
        throw std::invalid_argument(
            "Reference tetrahedron coordinates must be finite."
        );
    }
    TetrahedronP1BasisEvaluation3D result;
    result.values = {1.0 - r - s - t, r, s, t};
    result.reference_gradients = {{
        {{-1.0, -1.0, -1.0}},
        {{ 1.0,  0.0,  0.0}},
        {{ 0.0,  1.0,  0.0}},
        {{ 0.0,  0.0,  1.0}}
    }};
    return result;
}

TetrahedronP1PhysicalGradients3D
evaluate_tetrahedron_p1_physical_gradients_3d(
    const TetrahedralMesh3D& mesh,
    const std::size_t element
) {
    mesh.validate();
    const Tetrahedron3D& tetrahedron = mesh.elements.at(element);
    const Point3D& p0 = mesh.nodes[tetrahedron[0]];
    const Vector3 e1 = subtract(mesh.nodes[tetrahedron[1]], p0);
    const Vector3 e2 = subtract(mesh.nodes[tetrahedron[2]], p0);
    const Vector3 e3 = subtract(mesh.nodes[tetrahedron[3]], p0);
    const Vector3 c23 = cross(e2, e3);
    const Vector3 c31 = cross(e3, e1);
    const Vector3 c12 = cross(e1, e2);
    const double determinant = dot(e1, c23);
    if (!std::isfinite(determinant) || determinant == 0.0) {
        throw std::invalid_argument(
            "Cannot evaluate physical P1 gradients on a degenerate tetrahedron."
        );
    }

    const auto reference = evaluate_tetrahedron_p1_basis_3d(0.0, 0.0, 0.0);
    TetrahedronP1PhysicalGradients3D result;
    result.signed_jacobian_determinant = determinant;
    for (std::size_t local = 0; local < 4U; ++local) {
        const double dr = reference.reference_gradients[local][0];
        const double ds = reference.reference_gradients[local][1];
        const double dt = reference.reference_gradients[local][2];
        result.gradients[local] = {{
            (dr * c23.x + ds * c31.x + dt * c12.x) / determinant,
            (dr * c23.y + ds * c31.y + dt * c12.y) / determinant,
            (dr * c23.z + ds * c31.z + dt * c12.z) / determinant
        }};
    }
    return result;
}

} // namespace nova::fem
