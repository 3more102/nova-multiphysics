#include "p1_integration_3d.hpp"

#include <limits>

namespace nova::fem::detail {
namespace {

// K = D R D, with D_ii = sqrt(K_ii). Factor the unit-diagonal R so the
// pivot test is independent of material units and diagonal contrast.
void transform_tensor_gradients(const SymmetricTensor3D& k, P1Gradients3D& gradients) {
    if (!(k.xx > 0 && k.yy > 0 && k.zz > 0) ||
        !std::isfinite(k.xx) || !std::isfinite(k.yy) || !std::isfinite(k.zz)) {
        throw std::invalid_argument("Conductivity tensor diagonals must be finite and positive.");
    }
    const double dx = std::sqrt(k.xx), dy = std::sqrt(k.yy), dz = std::sqrt(k.zz);
    const auto correlation = [](double entry, double a, double b) {
        // Divide by the smaller root first to preserve weak cross terms on
        // axes with very different conductivities.
        return (entry / std::min(a, b)) / std::max(a, b);
    };
    const double xy = correlation(k.xy, dx, dy);
    const double xz = correlation(k.xz, dx, dz);
    const double yz = correlation(k.yz, dy, dz);
    if (!(std::abs(xy) < 1 && std::abs(xz) < 1 && std::abs(yz) < 1)) {
        throw std::invalid_argument("Conductivity tensor has invalid cross terms.");
    }
    constexpr double pivot_floor = 32 * std::numeric_limits<double>::epsilon();
    const double pivot_y = std::fma(-xy, xy, 1.0);
    if (pivot_y <= pivot_floor) {
        throw std::invalid_argument("Conductivity tensor has an unresolved second pivot.");
    }
    const double ly = std::sqrt(pivot_y);
    const double lzy = std::fma(-xy, xz, yz) / ly;
    const double pivot_z = std::fma(-lzy, lzy, std::fma(-xz, xz, 1.0));
    if (!(pivot_z > pivot_floor)) {
        throw std::invalid_argument("Conductivity tensor must be numerically positive definite.");
    }
    const double lz = std::sqrt(pivot_z);
    for (auto& g : gradients) {
        const double x = dx * g[0], y = dy * g[1], z = dz * g[2];
        g = {x + xy * y + xz * z, ly * y + lzy * z, lz * z};
    }
}

} // namespace

void ConductivityField3D::validate() const {
    if (static_cast<bool>(scalar) == static_cast<bool>(tensor)) {
        throw std::invalid_argument("Supply exactly one scalar or tensor 3D conductivity field.");
    }
}

void add_conductivity_sample(P1Matrix3D& stiffness, const ConductivityField3D& material,
                             const Point3D& point, const P1Gradients3D& gradients,
                             double sqrt_weight) {
    if (material.scalar) {
        const double k = material.scalar(point.x, point.y, point.z);
        if (!std::isfinite(k) || k <= 0) {
            throw std::invalid_argument("Scalar 3D conductivity must be finite and positive.");
        }
        sqrt_weight *= std::sqrt(k);
    }
    auto weighted = gradients;
    for (auto& g : weighted) {
        for (auto& component : g) component *= sqrt_weight;
    }
    if (material.tensor) {
        transform_tensor_gradients(material.tensor(point.x, point.y, point.z), weighted);
    }
    // Assemble a Gram matrix from sqrt(weight) L^T grad(N), K = L L^T.
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const auto& a = weighted[i];
            const auto& b = weighted[j];
            stiffness[i][j] += a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
            if (!std::isfinite(stiffness[i][j])) {
                throw std::overflow_error("3D element stiffness exceeds the supported numeric range.");
            }
        }
    }
}

} // namespace nova::fem::detail
