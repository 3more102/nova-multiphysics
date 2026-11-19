#include "nova/fem/error_norms_3d.hpp"

#include "nova/fem/tetrahedron_h1.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nova::fem {
namespace {

struct GaussPoint1D {
    double coordinate;
    double weight;
};

// Five-point Gauss-Legendre quadrature on [0,1], degree-nine exact.
constexpr std::array<GaussPoint1D, 5> gauss_legendre_5{{
    {0.046910077030668, 0.118463442528095},
    {0.230765344947158, 0.239314335249683},
    {0.5,               0.284444444444444},
    {0.769234655052842, 0.239314335249683},
    {0.953089922969332, 0.118463442528095}
}};

[[nodiscard]] double checked_norm(const long double square) {
    if (!(square >= 0.0L) || !std::isfinite(square)) {
        throw std::overflow_error("3D error-norm integration is non-finite.");
    }
    const long double magnitude = std::sqrt(square);
    if (!std::isfinite(magnitude) ||
        magnitude > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("3D error norm exceeds the finite range.");
    }
    return static_cast<double>(magnitude);
}

} // namespace

ErrorNorms3D compute_error_norms_3d(
    const TetrahedralMesh3D& mesh,
    const std::span<const double> nodal_values,
    const ExactScalarField3D& exact_solution
) {
    mesh.validate();
    if (nodal_values.size() != mesh.nodes.size()) {
        throw std::invalid_argument("3D error-norm field size must match node count.");
    }
    if (!exact_solution.value || !exact_solution.derivative_x ||
        !exact_solution.derivative_y || !exact_solution.derivative_z) {
        throw std::invalid_argument("3D exact value and all three derivatives are required.");
    }
    for (const double value : nodal_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("3D error-norm nodal field must be finite.");
        }
    }

    long double l2_squared = 0.0L;
    long double h1_seminorm_squared = 0.0L;
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const auto& cell = mesh.elements[e];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, e);
        const double jacobian = std::abs(basis.signed_jacobian_determinant);
        const auto& p0 = mesh.nodes[cell[0]];
        const auto& p1 = mesh.nodes[cell[1]];
        const auto& p2 = mesh.nodes[cell[2]];
        const auto& p3 = mesh.nodes[cell[3]];
        std::array<double, 4> values{};
        std::array<double, 3> numerical_gradient{};
        for (std::size_t a = 0; a < 4; ++a) {
            values[a] = nodal_values[cell[a]];
            for (std::size_t d = 0; d < 3; ++d) {
                numerical_gradient[d] += values[a] * basis.gradients[a][d];
            }
        }
        for (const auto& q_u : gauss_legendre_5) {
            const double u = q_u.coordinate;
            const double one_minus_u = 1.0 - u;
            for (const auto& q_v : gauss_legendre_5) {
                const double v = q_v.coordinate;
                const double one_minus_v = 1.0 - v;
                for (const auto& q_w : gauss_legendre_5) {
                    // Duffy map [0,1]^3 -> reference tetrahedron:
                    // r=u, s=(1-u)v, t=(1-u)(1-v)w.
                    const double r = u;
                    const double s = one_minus_u * v;
                    const double t = one_minus_u * one_minus_v * q_w.coordinate;
                    const std::array<double, 4> shape{
                        1.0 - r - s - t, r, s, t
                    };
                    const double x = p0.x + r * (p1.x - p0.x) +
                        s * (p2.x - p0.x) + t * (p3.x - p0.x);
                    const double y = p0.y + r * (p1.y - p0.y) +
                        s * (p2.y - p0.y) + t * (p3.y - p0.y);
                    const double z = p0.z + r * (p1.z - p0.z) +
                        s * (p2.z - p0.z) + t * (p3.z - p0.z);
                    double numerical_value = 0.0;
                    for (std::size_t a = 0; a < 4; ++a) {
                        numerical_value += values[a] * shape[a];
                    }
                    const double exact_value = exact_solution.value(x, y, z);
                    const std::array<double, 3> exact_gradient{
                        exact_solution.derivative_x(x, y, z),
                        exact_solution.derivative_y(x, y, z),
                        exact_solution.derivative_z(x, y, z)
                    };
                    if (!std::isfinite(exact_value) ||
                        !std::isfinite(exact_gradient[0]) ||
                        !std::isfinite(exact_gradient[1]) ||
                        !std::isfinite(exact_gradient[2])) {
                        throw std::invalid_argument("3D exact field returned non-finite data.");
                    }

                    const long double error =
                        static_cast<long double>(numerical_value) - exact_value;
                    long double gradient_error_squared = 0.0L;
                    for (std::size_t d = 0; d < 3; ++d) {
                        const long double component =
                            static_cast<long double>(numerical_gradient[d]) -
                            exact_gradient[d];
                        gradient_error_squared += component * component;
                    }
                    const long double weight =
                        static_cast<long double>(jacobian) *
                        q_u.weight * q_v.weight * q_w.weight *
                        one_minus_u * one_minus_u * one_minus_v;
                    l2_squared += weight * error * error;
                    h1_seminorm_squared += weight * gradient_error_squared;
                }
            }
        }
        if (!std::isfinite(l2_squared) || !std::isfinite(h1_seminorm_squared)) {
            throw std::overflow_error("3D error-norm accumulation overflowed.");
        }
    }
    const double l2 = checked_norm(l2_squared);
    const double semi = checked_norm(h1_seminorm_squared);
    return {l2, semi, checked_norm(l2_squared + h1_seminorm_squared)};
}

} // namespace nova::fem
