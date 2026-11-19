#include "p1_integration_3d.hpp"

#include "nova/fem/diffusion_3d.hpp"

#include <limits>

namespace nova::fem::detail {
namespace {

// K = D R D, with D_ii = sqrt(K_ii). Factor the unit-diagonal R so the
// pivot test is independent of material units and diagonal contrast.
struct TensorFactor {
    double dx, dy, dz;
    double xy, xz, ly, lzy, lz;

    void transform(P1Gradients3D& gradients) const {
        for (auto& g : gradients) {
            const double x = dx * g[0], y = dy * g[1], z = dz * g[2];
            g = {x + xy * y + xz * z, ly * y + lzy * z, lz * z};
        }
    }
};

[[nodiscard]] TensorFactor factor_tensor(const SymmetricTensor3D& k) {
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
    return {dx, dy, dz, xy, xz, ly, lzy, std::sqrt(pivot_z)};
}

void check_scalar(double k) {
    if (!std::isfinite(k) || k <= 0.0) {
        throw std::invalid_argument("Scalar 3D conductivity must be finite and positive.");
    }
}

} // namespace

void ConductivityField3D::validate() const {
    const int choices = (scalar ? 1 : 0) + (tensor ? 1 : 0) + (constant ? 1 : 0);
    if (choices != 1) {
        throw std::invalid_argument("Select exactly one 3D conductivity representation.");
    }
}

double ConductivityField3D::scalar_at(const Point3D& point) const {
    if (constant) return *constant; // Checked once when the cell table is bound.
    const double k = scalar(point.x, point.y, point.z);
    check_scalar(k);
    return k;
}

Point3D ConductivityField3D::diffusive_flux(const Point3D& point, const Point3D& g) const {
    Point3D flux;
    if (tensor) {
        const auto k = tensor(point.x, point.y, point.z);
        static_cast<void>(factor_tensor(k));
        flux = {-(k.xx*g.x + k.xy*g.y + k.xz*g.z),
                -(k.xy*g.x + k.yy*g.y + k.yz*g.z),
                -(k.xz*g.x + k.yz*g.y + k.zz*g.z)};
    } else {
        const double k = scalar_at(point);
        flux = {-k*g.x, -k*g.y, -k*g.z};
    }
    if (!std::isfinite(flux.x) || !std::isfinite(flux.y) || !std::isfinite(flux.z)) {
        throw std::overflow_error("3D diffusive flux exceeds the finite range.");
    }
    return flux;
}

DiffusionMaterials3D::DiffusionMaterials3D(const ScalarDiffusion3DProblem& problem)
    : problem_(problem) {
    if (problem.cell_diffusivity.empty()) {
        ConductivityField3D{problem.diffusivity, problem.diffusivity_tensor}.validate();
        return;
    }
    if (problem.diffusivity || problem.diffusivity_tensor) {
        throw std::invalid_argument("Cell conductivity cannot be combined with a spatial field.");
    }
    if (problem.cell_diffusivity.size() != problem.mesh.elements.size()) {
        throw std::invalid_argument("Cell conductivity requires one value per tetrahedron.");
    }
    for (const double k : problem.cell_diffusivity) check_scalar(k);
}

ConductivityField3D DiffusionMaterials3D::on_cell(std::size_t element) const {
    ConductivityField3D material{problem_.diffusivity, problem_.diffusivity_tensor};
    if (!problem_.cell_diffusivity.empty()) material.constant = problem_.cell_diffusivity[element];
    return material;
}

void add_conductivity_sample(P1Matrix3D& stiffness, const ConductivityField3D& material,
                             const Point3D& point, const P1Gradients3D& gradients,
                             double sqrt_weight) {
    if (!material.tensor) {
        sqrt_weight *= std::sqrt(material.scalar_at(point));
    }
    auto weighted = gradients;
    for (auto& g : weighted) {
        for (auto& component : g) component *= sqrt_weight;
    }
    if (material.tensor) {
        factor_tensor(material.tensor(point.x, point.y, point.z)).transform(weighted);
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
