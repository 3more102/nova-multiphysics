#include "nova/fem/p1_flux_3d.hpp"

#include "p1_integration_3d.hpp"

#include <cmath>
#include <stdexcept>

namespace nova::fem {
namespace {

[[nodiscard]] bool finite(const Point3D& value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
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
    const detail::DiffusionMaterials3D materials(problem);

    std::vector<P1TetrahedralFlux3D> output;
    output.reserve(problem.mesh.elements.size());
    for (std::size_t index = 0; index < problem.mesh.elements.size(); ++index) {
        const auto& tetra = problem.mesh.elements[index];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(
            problem.mesh, index
        );
        Point3D gradient{};
        // Partition of unity lets us subtract a reference temperature before
        // differentiating. A uniform field then has exactly zero flux.
        for (std::size_t local = 1; local < tetra.size(); ++local) {
            const double node = nodal_values[tetra[local]];
            const double origin = nodal_values[tetra[0]];
            double difference = node - origin;
            double scale = 1.0;
            if (!std::isfinite(difference)) {
                // Opposite extreme values can still give a finite gradient on a long edge.
                difference = 0.5 * node - 0.5 * origin;
                scale = 2.0;
            }
            gradient.x += (difference * basis.gradients[local][0]) * scale;
            gradient.y += (difference * basis.gradients[local][1]) * scale;
            gradient.z += (difference * basis.gradients[local][2]) * scale;
        }
        if (!finite(gradient)) {
            throw std::overflow_error("3D P1 scalar gradient exceeds the finite range.");
        }

        const auto material = materials.on_cell(index);
        if (material.constant) {
            output.push_back({gradient, material.diffusive_flux(problem.mesh.nodes[tetra[0]], gradient)});
            continue;
        }
        Point3D average_flux{};
        // Four positive degree-two samples with weights summing to 1/6.
        // Multiplying by six returns the cell-volume average of q = -K grad(u).
        for (const auto& quadrature : tetrahedron_degree_two_quadrature) {
            const Point3D at = detail::simplex_point(
                problem.mesh, tetra, quadrature.shape
            );
            const Point3D sample = material.diffusive_flux(at, gradient);
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
