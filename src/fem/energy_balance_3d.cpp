#include "nova/fem/energy_balance_3d.hpp"

#include "p1_integration_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <vector>

namespace nova::fem {
namespace {

void prescribe(std::vector<std::optional<double>>& prescribed,
               std::size_t node, double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Energy balance requires finite Dirichlet values.");
    }
    if (prescribed[node]) {
        const double original = *prescribed[node];
        const double scale = std::max({1.0, std::abs(original), std::abs(value)});
        if (!(std::abs(original - value) <=
              64.0 * std::numeric_limits<double>::epsilon() * scale)) {
            throw std::invalid_argument("Energy balance has conflicting Dirichlet values.");
        }
        return;
    }
    prescribed[node] = value;
}

double checked(long double amount, const char* message) {
    if (!std::isfinite(amount) ||
        std::abs(amount) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<double>(amount);
}

} // namespace

DiffusionEnergyBalance3D evaluate_diffusion_energy_balance_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> nodal_values
) {
    detail::ExteriorFaceAssignments3D exterior(problem.mesh);
    const auto& mesh = problem.mesh;
    const std::size_t n = mesh.nodes.size();
    if (nodal_values.size() != n) {
        throw std::invalid_argument("Energy balance requires one finite temperature per mesh vertex.");
    }
    for (double value : nodal_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Energy balance nodal temperatures must be finite.");
        }
    }
    if (!problem.cell_source.empty()) {
        if (problem.cell_source.size() != mesh.elements.size()) {
            throw std::invalid_argument("Energy balance cell source must cover every tetrahedron.");
        }
        for (double value : problem.cell_source) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("Energy balance cell source must be finite.");
            }
        }
    }
    const detail::DiffusionMaterials3D materials(problem);
    std::vector<std::optional<double>> prescribed(n);
    std::set<std::size_t> explicit_nodes;
    for (const auto& constraint : problem.dirichlet_constraints) {
        if (constraint.node >= n) {
            throw std::out_of_range("Energy balance Dirichlet node does not exist.");
        }
        if (!explicit_nodes.insert(constraint.node).second) {
            throw std::invalid_argument("Energy balance duplicate Dirichlet vertex.");
        }
        prescribe(prescribed, constraint.node, constraint.value);
    }
    for (const auto& condition : problem.dirichlet_faces) {
        exterior.claim(condition.face);
        if (!condition.value) {
            throw std::invalid_argument("Energy balance Dirichlet face has no callback.");
        }
        for (std::size_t node : condition.face) {
            const auto& p = mesh.nodes[node];
            prescribe(prescribed, node, condition.value(p.x, p.y, p.z));
        }
    }
    for (const auto& condition : problem.neumann_faces) {
        exterior.claim(condition.face);
        if (!condition.outward_flux) {
            throw std::invalid_argument("Energy balance Neumann face has no callback.");
        }
    }
    for (const auto& condition : problem.robin_faces) {
        exterior.claim(condition.face);
        if (!condition.transfer_coefficient || !condition.ambient_value) {
            throw std::invalid_argument("Energy balance Robin face has missing callbacks.");
        }
    }
    for (std::size_t node = 0; node < n; ++node) {
        if (!prescribed[node]) continue;
        const double value = nodal_values[node];
        const double fixed = *prescribed[node];
        const double scale = std::max({1.0, std::abs(fixed), std::abs(value)});
        if (!(std::abs(value - fixed) <=
              64.0 * std::numeric_limits<double>::epsilon() * scale)) {
            throw std::invalid_argument("Energy balance field disagrees with fixed Dirichlet temperature.");
        }
    }

    std::vector<long double> residual(n, 0.0L);
    long double generation = 0.0L;
    long double outward_neumann = 0.0L;
    long double outward_robin = 0.0L;

    // Reuse precisely the solver's element integration rule and material
    // selection. Accumulate residual Au - b in original, unconstrained DOFs.
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const auto& tetra = mesh.elements[element];
        const auto physical = evaluate_tetrahedron_p1_physical_gradients_3d(mesh, element);
        const double local_source = problem.cell_source.empty()
            ? 0.0 : problem.cell_source[element];
        const auto local_matrix = detail::integrate_p1_volume_3d(
            mesh, tetra, physical, materials.on_cell(element),
            [&](const Point3D& point, const auto& shape, double weight) {
                const double spatial = problem.source
                    ? problem.source(point.x, point.y, point.z) : 0.0;
                if (!std::isfinite(spatial)) {
                    throw std::invalid_argument("Energy balance source callback is nonfinite.");
                }
                const double density = spatial + local_source;
                if (!std::isfinite(density)) {
                    throw std::overflow_error("Energy balance combined source density overflowed.");
                }
                generation += static_cast<long double>(density) * weight;
                for (std::size_t i = 0; i < 4; ++i) {
                    const double load = density * shape[i] * weight;
                    if (!std::isfinite(load)) {
                        throw std::overflow_error("Energy balance element load overflowed.");
                    }
                    residual[tetra[i]] -= load;
                }
            });
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = 0; j < 4; ++j) {
                residual[tetra[i]] +=
                    static_cast<long double>(local_matrix[i][j]) * nodal_values[tetra[j]];
            }
        }
    }

    for (const auto& condition : problem.neumann_faces) {
        const auto& face = condition.face;
        const double jacobian = detail::face_jacobian(mesh, face);
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto point = detail::simplex_point(mesh, face, q.shape);
            const double outward = condition.outward_flux(point.x, point.y, point.z);
            if (!std::isfinite(outward)) {
                throw std::invalid_argument("Energy balance Neumann flux is nonfinite.");
            }
            const double weight = jacobian * q.weight;
            outward_neumann += static_cast<long double>(outward) * weight;
            for (std::size_t i = 0; i < 3; ++i) {
                const double load = outward * q.shape[i] * weight;
                if (!std::isfinite(load)) {
                    throw std::overflow_error("Energy balance Neumann load overflowed.");
                }
                residual[face[i]] += load;
            }
        }
    }

    for (const auto& condition : problem.robin_faces) {
        const auto& face = condition.face;
        const double jacobian = detail::face_jacobian(mesh, face);
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto point = detail::simplex_point(mesh, face, q.shape);
            const double h = condition.transfer_coefficient(point.x, point.y, point.z);
            const double ambient = condition.ambient_value(point.x, point.y, point.z);
            if (!std::isfinite(h) || h < 0.0 || !std::isfinite(ambient)) {
                throw std::invalid_argument("Energy balance Robin coefficients are invalid.");
            }
            const double ambient_load = h * ambient;
            if (!std::isfinite(ambient_load)) {
                throw std::overflow_error("Energy balance Robin ambient load overflowed.");
            }
            const double weight = jacobian * q.weight;
            long double temperature = 0.0L;
            for (std::size_t i = 0; i < 3; ++i) {
                temperature += static_cast<long double>(nodal_values[face[i]]) * q.shape[i];
            }
            outward_robin += static_cast<long double>(h) *
                             (temperature - ambient) * weight;
            for (std::size_t i = 0; i < 3; ++i) {
                const double load = ambient_load * q.shape[i] * weight;
                if (!std::isfinite(load)) {
                    throw std::overflow_error("Energy balance Robin load overflowed.");
                }
                residual[face[i]] -= load;
                for (std::size_t j = 0; j < 3; ++j) {
                    const double matrix = h * q.shape[i] * q.shape[j] * weight;
                    if (!std::isfinite(matrix)) {
                        throw std::overflow_error("Energy balance Robin matrix overflowed.");
                    }
                    residual[face[i]] += static_cast<long double>(matrix) * nodal_values[face[j]];
                }
            }
        }
    }

    long double constrained_residual = 0.0L;
    long double max_free_residual = 0.0L;
    std::size_t fixed_count = 0;
    for (std::size_t node = 0; node < n; ++node) {
        if (!std::isfinite(residual[node])) {
            throw std::overflow_error("Energy balance weak residual is nonfinite.");
        }
        if (prescribed[node]) {
            ++fixed_count;
            constrained_residual += residual[node];
        } else {
            max_free_residual = std::max(max_free_residual, std::abs(residual[node]));
        }
    }
    const long double dirichlet_outward = -constrained_residual;
    const long double defect = generation - outward_neumann -
                               outward_robin - dirichlet_outward;
    return {
        checked(generation, "Generated heat exceeds finite range."),
        checked(outward_neumann, "Neumann heat exceeds finite range."),
        checked(outward_robin, "Robin heat exceeds finite range."),
        checked(dirichlet_outward, "Constraint reaction exceeds finite range."),
        checked(defect, "Energy balance defect exceeds finite range."),
        checked(max_free_residual, "Free equation residual exceeds finite range."),
        fixed_count, n - fixed_count
    };
}

} // namespace nova::fem
