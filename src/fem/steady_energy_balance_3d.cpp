#include "nova/fem/steady_energy_balance_3d.hpp"

#include "p1_integration_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nova::fem {
namespace {

[[nodiscard]] double finite_rate(const long double value) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("Steady energy-balance rate is outside the finite double range.");
    }
    return static_cast<double>(value);
}

} // namespace

SteadyEnergyBalance3D evaluate_steady_energy_balance_3d(
    const ScalarDiffusion3DProblem& problem,
    const std::span<const double> temperatures
) {
    detail::ExteriorFaceAssignments3D exterior(problem.mesh);
    const std::size_t node_count = problem.mesh.nodes.size();
    if (temperatures.size() != node_count) {
        throw std::invalid_argument("Energy balance requires one temperature per mesh vertex.");
    }
    for (const double value : temperatures) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Energy balance requires finite nodal temperatures.");
        }
    }
    if (!problem.source) {
        throw std::invalid_argument("Energy balance requires a source field.");
    }
    if (!problem.cell_source.empty() &&
        problem.cell_source.size() != problem.mesh.elements.size()) {
        throw std::invalid_argument("Energy balance requires one cell source per tetrahedron.");
    }
    for (const double value : problem.cell_source) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Energy balance requires finite cell source densities.");
        }
    }
    const detail::DiffusionMaterials3D materials(problem);

    std::vector<unsigned char> is_dirichlet(node_count, 0U);
    for (const auto& condition : problem.dirichlet_constraints) {
        if (condition.node >= node_count || !std::isfinite(condition.value)) {
            throw std::invalid_argument("Invalid Dirichlet vertex in energy balance.");
        }
        is_dirichlet[condition.node] = 1U;
    }
    for (const auto& condition : problem.dirichlet_faces) {
        exterior.claim(condition.face);
        if (!condition.value) {
            throw std::invalid_argument("Missing Dirichlet face field in energy balance.");
        }
        for (const std::size_t node : condition.face) is_dirichlet[node] = 1U;
    }
    for (const auto& condition : problem.neumann_faces) {
        exterior.claim(condition.face);
        if (!condition.outward_flux) {
            throw std::invalid_argument("Missing Neumann face field in energy balance.");
        }
    }
    for (const auto& condition : problem.robin_faces) {
        exterior.claim(condition.face);
        if (!condition.transfer_coefficient || !condition.ambient_value) {
            throw std::invalid_argument("Missing Robin face field in energy balance.");
        }
    }

    long double source_rate = 0.0L;
    long double reaction = 0.0L;
    long double neumann_rate = 0.0L;
    long double robin_rate = 0.0L;

    for (std::size_t e = 0; e < problem.mesh.elements.size(); ++e) {
        const auto& cell = problem.mesh.elements[e];
        const auto basis = evaluate_tetrahedron_p1_physical_gradients_3d(problem.mesh, e);
        const double cell_source = problem.cell_source.empty() ? 0.0 : problem.cell_source[e];
        const auto stiffness = detail::integrate_p1_volume_3d(
            problem.mesh, cell, basis, materials.on_cell(e),
            [&](const Point3D& point, const auto& shape, const double weight) {
                const double spatial = problem.source(point.x, point.y, point.z);
                if (!std::isfinite(spatial)) {
                    throw std::invalid_argument("Nonfinite spatial source in energy balance.");
                }
                const double density = spatial + cell_source;
                if (!std::isfinite(density)) {
                    throw std::overflow_error("Combined source exceeds finite range.");
                }
                source_rate += static_cast<long double>(density) * weight;
                for (std::size_t i = 0; i < 4; ++i) {
                    if (is_dirichlet[cell[i]]) {
                        // Same order of floating-point operations as the FEM load.
                        reaction -= static_cast<long double>(
                            density * shape[i] * weight);
                    }
                }
            });
        for (std::size_t i = 0; i < 4; ++i) {
            if (!is_dirichlet[cell[i]]) continue;
            for (std::size_t j = 0; j < 4; ++j) {
                reaction += static_cast<long double>(stiffness[i][j]) *
                            static_cast<long double>(temperatures[cell[j]]);
            }
        }
    }

    for (const auto& condition : problem.neumann_faces) {
        const double jacobian = detail::face_jacobian(problem.mesh, condition.face);
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto point = detail::simplex_point(problem.mesh, condition.face, q.shape);
            const double flux = condition.outward_flux(point.x, point.y, point.z);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument("Nonfinite Neumann flux in energy balance.");
            }
            const double weight = jacobian * q.weight;
            neumann_rate += static_cast<long double>(flux) * weight;
            for (std::size_t i = 0; i < 3; ++i) {
                if (is_dirichlet[condition.face[i]]) {
                    reaction += static_cast<long double>(
                        flux * q.shape[i] * weight);
                }
            }
        }
    }

    for (const auto& condition : problem.robin_faces) {
        const double jacobian = detail::face_jacobian(problem.mesh, condition.face);
        for (const auto& q : triangle_face_degree_two_quadrature) {
            const auto point = detail::simplex_point(problem.mesh, condition.face, q.shape);
            const double h = condition.transfer_coefficient(point.x, point.y, point.z);
            const double ambient = condition.ambient_value(point.x, point.y, point.z);
            if (!std::isfinite(h) || h < 0.0 || !std::isfinite(ambient)) {
                throw std::invalid_argument("Invalid Robin data in energy balance.");
            }
            long double temperature = 0.0L;
            for (std::size_t i = 0; i < 3; ++i) {
                temperature += static_cast<long double>(q.shape[i]) *
                    static_cast<long double>(temperatures[condition.face[i]]);
            }
            const long double outward = static_cast<long double>(h) *
                (temperature - static_cast<long double>(ambient));
            const double weight = jacobian * q.weight;
            robin_rate += outward * weight;
            for (std::size_t i = 0; i < 3; ++i) {
                if (is_dirichlet[condition.face[i]]) {
                    reaction += outward * q.shape[i] * weight;
                }
            }
        }
    }

    const long double dirichlet_rate = -reaction;
    const long double outward_rate = dirichlet_rate + neumann_rate + robin_rate;
    return {
        finite_rate(source_rate),
        finite_rate(neumann_rate),
        finite_rate(robin_rate),
        finite_rate(dirichlet_rate),
        finite_rate(outward_rate),
        finite_rate(outward_rate - source_rate)
    };
}

} // namespace nova::fem
