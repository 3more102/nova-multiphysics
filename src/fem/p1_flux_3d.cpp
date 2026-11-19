#include "nova/fem/p1_flux_3d.hpp"

#include "p1_integration_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace nova::fem {
namespace {

[[nodiscard]] bool finite(const Point3D& value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
}

// Match the normalized SPD acceptance criteria of the tensor stiffness
// assembly, so postprocessing cannot silently accept an invalid material.
void validate_tensor(const SymmetricTensor3D& k) {
    if (!(k.xx > 0.0 && k.yy > 0.0 && k.zz > 0.0) ||
        !std::isfinite(k.xx) || !std::isfinite(k.yy) ||
        !std::isfinite(k.zz) || !std::isfinite(k.xy) ||
        !std::isfinite(k.xz) || !std::isfinite(k.yz)) {
        throw std::invalid_argument("3D flux requires a finite tensor with positive diagonal.");
    }
    const double dx = std::sqrt(k.xx);
    const double dy = std::sqrt(k.yy);
    const double dz = std::sqrt(k.zz);
    const auto normalized = [](double value, double a, double b) {
        return (value / std::min(a, b)) / std::max(a, b);
    };
    const double xy = normalized(k.xy, dx, dy);
    const double xz = normalized(k.xz, dx, dz);
    const double yz = normalized(k.yz, dy, dz);
    if (!(std::abs(xy) < 1.0 && std::abs(xz) < 1.0 && std::abs(yz) < 1.0)) {
        throw std::invalid_argument("3D flux tensor has invalid cross terms.");
    }
    constexpr double floor = 32.0 * std::numeric_limits<double>::epsilon();
    const double pivot_y = std::fma(-xy, xy, 1.0);
    if (!(pivot_y > floor)) {
        throw std::invalid_argument("3D flux tensor has an unresolved second pivot.");
    }
    const double lzy = std::fma(-xy, xz, yz) / std::sqrt(pivot_y);
    const double pivot_z = std::fma(-lzy, lzy, std::fma(-xz, xz, 1.0));
    if (!(pivot_z > floor)) {
        throw std::invalid_argument("3D flux tensor must be numerically positive definite.");
    }
}

[[nodiscard]] Point3D negative_material_gradient(
    const ScalarDiffusion3DProblem& problem,
    const Point3D& at,
    const Point3D& gradient
) {
    Point3D flux{};
    if (problem.diffusivity) {
        const double k = problem.diffusivity(at.x, at.y, at.z);
        if (!(k > 0.0) || !std::isfinite(k)) {
            throw std::invalid_argument("3D flux scalar diffusivity must be finite and positive.");
        }
        flux = {-k * gradient.x, -k * gradient.y, -k * gradient.z};
    } else {
        const SymmetricTensor3D k = problem.diffusivity_tensor(at.x, at.y, at.z);
        validate_tensor(k);
        flux = {
            -(k.xx * gradient.x + k.xy * gradient.y + k.xz * gradient.z),
            -(k.xy * gradient.x + k.yy * gradient.y + k.yz * gradient.z),
            -(k.xz * gradient.x + k.yz * gradient.y + k.zz * gradient.z)
        };
    }
    if (!finite(flux)) {
        throw std::overflow_error("3D diffusive flux exceeds the finite range.");
    }
    return flux;
}

} // namespace

std::vector<P1TetrahedralFlux3D> evaluate_p1_tetrahedral_flux_3d(
    const ScalarDiffusion3DProblem& problem,
    const std::span<const double> nodal_values
) {
    problem.mesh.validate();
    if (nodal_values.size() != problem.mesh.nodes.size()) {
        throw std::invalid_argument("3D flux requires one scalar value per mesh vertex.");
    }
    for (const double value : nodal_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("3D flux nodal values must be finite.");
        }
    }
    detail::ConductivityField3D{problem.diffusivity, problem.diffusivity_tensor}.validate();

    std::vector<P1TetrahedralFlux3D> output;
    output.reserve(problem.mesh.elements.size());
    for (std::size_t index = 0; index < problem.mesh.elements.size(); ++index) {
        const auto& tetra = problem.mesh.elements[index];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(
            problem.mesh, index
        );
        Point3D gradient{};
        for (std::size_t local = 0; local < tetra.size(); ++local) {
            const double value = nodal_values[tetra[local]];
            gradient.x += value * basis.gradients[local][0];
            gradient.y += value * basis.gradients[local][1];
            gradient.z += value * basis.gradients[local][2];
        }
        if (!finite(gradient)) {
            throw std::overflow_error("3D P1 scalar gradient exceeds the finite range.");
        }

        Point3D average_flux{};
        // Four positive degree-two samples with weights summing to 1/6.
        // Multiplying by six returns the cell-volume average of q = -K grad(u).
        for (const auto& quadrature : tetrahedron_degree_two_quadrature) {
            const Point3D at = detail::simplex_point(
                problem.mesh, tetra, quadrature.shape
            );
            const Point3D sample = negative_material_gradient(problem, at, gradient);
            const double weight = 6.0 * quadrature.weight;
            average_flux.x += weight * sample.x;
            average_flux.y += weight * sample.y;
            average_flux.z += weight * sample.z;
        }
        if (!finite(average_flux)) {
            throw std::overflow_error("3D P1 element-average flux is non-finite.");
        }
        output.push_back({gradient, average_flux});
    }
    return output;
}

} // namespace nova::fem
