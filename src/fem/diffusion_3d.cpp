#include "nova/fem/diffusion_3d.hpp"

#include "p1_integration_3d.hpp"
#include "diffusion_nullspace.hpp"
#include "constrained_diffusion.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

struct RobinEnergySample {
    BoundaryFace3D face{};
    std::array<double, 3> shape{};
    double transfer{};
    double ambient{};
    double weight{};
};

[[nodiscard]] double finite_energy_rate(const long double rate) {
    if (!std::isfinite(rate) ||
        std::abs(rate) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error("3D energy balance is outside finite double range.");
    }
    return static_cast<double>(rate);
}

void prescribe(
    std::vector<std::optional<double>>& prescribed,
    const std::size_t node,
    const double value
) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("Dirichlet values must be finite.");
    }
    if (prescribed[node].has_value()) {
        const double previous = *prescribed[node];
        const double scale = std::max({
            1.0, std::abs(previous), std::abs(value)
        });
        if (std::abs(previous - value) >
            64.0 * std::numeric_limits<double>::epsilon() * scale) {
            throw std::invalid_argument(
                "Conflicting Dirichlet values at a shared 3D vertex."
            );
        }
        return;
    }
    prescribed[node] = value;
}

} // namespace

ScalarDiffusion3DResult solve_scalar_diffusion_3d(
    const ScalarDiffusion3DProblem& problem,
    const linalg::SolverOptions& solver_options
) {
    detail::ExteriorFaceAssignments3D exterior(problem.mesh);
    const std::size_t node_count = problem.mesh.nodes.size();
    detail::DiffusionNullspace nullspace(node_count, problem.mesh.elements.size(),
        [&](std::size_t e) -> std::span<const std::size_t> { return problem.mesh.elements[e]; });

    const detail::DiffusionMaterials3D materials(problem);
    if (!problem.source) throw std::invalid_argument("3D diffusion requires a source field.");
    if (!problem.cell_source.empty()) {
        if (problem.cell_source.size() != problem.mesh.elements.size()) {
            throw std::invalid_argument("Cell source requires one value per tetrahedron.");
        }
        for (const double value : problem.cell_source) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("Cell source densities must be finite.");
            }
        }
    }

    std::vector<std::optional<double>> prescribed(node_count);
    std::set<std::size_t> explicit_nodes;
    for (const DirichletConstraint3D& constraint :
         problem.dirichlet_constraints) {
        if (constraint.node >= node_count) {
            throw std::out_of_range(
                "3D Dirichlet constraint references a missing vertex."
            );
        }
        if (!explicit_nodes.insert(constraint.node).second) {
            throw std::invalid_argument(
                "Duplicate explicit 3D Dirichlet vertex constraint."
            );
        }
        prescribe(prescribed, constraint.node, constraint.value);
    }

    for (const DirichletBoundaryFace3D& condition :
         problem.dirichlet_faces) {
        exterior.claim(condition.face);
        if (!condition.value) {
            throw std::invalid_argument(
                "3D Dirichlet face requires a value function."
            );
        }
        for (const std::size_t node : condition.face) {
            const Point3D& point = problem.mesh.nodes[node];
            prescribe(
                prescribed, node,
                condition.value(point.x, point.y, point.z)
            );
        }
    }

    for (const NeumannBoundaryFace3D& condition : problem.neumann_faces) {
        exterior.claim(condition.face);
        if (!condition.outward_flux) {
            throw std::invalid_argument(
                "3D Neumann face requires an outward-flux function."
            );
        }
    }

    for (const RobinBoundaryFace3D& condition : problem.robin_faces) {
        exterior.claim(condition.face);
        if (!condition.transfer_coefficient || !condition.ambient_value) {
            throw std::invalid_argument(
                "3D Robin face requires transfer and ambient functions."
            );
        }
    }

    constexpr std::size_t local_dof_count = 4U;
    constexpr std::size_t volume_entries_per_element =
        local_dof_count * local_dof_count;
    constexpr std::size_t robin_matrix_entries_per_face = 3U * 3U * 3U;

    std::vector<linalg::Triplet> entries;
    const std::size_t capacity = entries.max_size();
    if (problem.mesh.elements.size() >
        capacity / volume_entries_per_element) {
        throw std::length_error("3D finite-element volume matrix is too large.");
    }
    const std::size_t volume_entry_count =
        problem.mesh.elements.size() * volume_entries_per_element;
    if (problem.robin_faces.size() >
        (capacity - volume_entry_count) / robin_matrix_entries_per_face) {
        throw std::length_error("3D Robin boundary matrix is too large.");
    }
    entries.reserve(
        volume_entry_count +
        problem.robin_faces.size() * robin_matrix_entries_per_face
    );

    std::vector<double> rhs(node_count, 0.0);
    long double generated_heat_rate = 0.0L;
    long double neumann_heat_rate = 0.0L;
    std::vector<RobinEnergySample> robin_energy_samples;

    for (std::size_t element = 0;
         element < problem.mesh.elements.size();
         ++element) {
        const Tetrahedron3D& cell = problem.mesh.elements[element];
        const auto physical = evaluate_tetrahedron_p1_physical_gradients_3d(
            problem.mesh, element
        );
        const double absolute_determinant =
            std::abs(physical.signed_jacobian_determinant);
        const long double nodal_volume_weight =
            static_cast<long double>(absolute_determinant) / 24.0L;
        for (const std::size_t node : cell) {
            nullspace.add_basis_integral(node, nodal_volume_weight);
        }

        std::array<double, 4> local_rhs{};
        const double cell_source = problem.cell_source.empty() ? 0.0 : problem.cell_source[element];
        const auto local_matrix = detail::integrate_p1_volume_3d(
            problem.mesh, cell, physical, materials.on_cell(element),
            [&](const Point3D& point, const auto& shape, double weight) {
                const double spatial_source = problem.source(point.x, point.y, point.z);
                if (!std::isfinite(spatial_source)) {
                    throw std::invalid_argument("3D source must be finite.");
                }
                const double source = spatial_source + cell_source;
                if (!std::isfinite(source)) {
                    throw std::overflow_error("Combined 3D source density exceeds the finite range.");
                }
                if (problem.compute_energy_balance) {
                    generated_heat_rate += static_cast<long double>(source) * weight;
                }
                for (std::size_t i = 0; i < local_dof_count; ++i) {
                    const double load = source * shape[i] * weight;
                    local_rhs[i] += load;
                    nullspace.record_load(cell[i], load);
                }
            });

        for (std::size_t i = 0; i < local_dof_count; ++i) {
            const std::size_t global_i = cell[i];
            rhs[global_i] += local_rhs[i];
            if (!std::isfinite(rhs[global_i])) {
                throw std::overflow_error(
                    "3D element load-vector accumulation overflowed."
                );
            }
            for (std::size_t j = 0; j < local_dof_count; ++j) {
                entries.push_back({
                    global_i, cell[j], local_matrix[i][j]
                });
            }
        }
    }

    const auto integrate_neumann_face = [&](
        const BoundaryFace3D& face,
        const ScalarField3D& outward_flux
    ) {
        const double area_scale = detail::face_jacobian(problem.mesh, face);
        for (const auto& quadrature : triangle_face_degree_two_quadrature) {
            const auto point = detail::simplex_point(problem.mesh, face, quadrature.shape);
            const double flux = outward_flux(point.x, point.y, point.z);
            if (!std::isfinite(flux)) {
                throw std::invalid_argument(
                    "3D Neumann outward flux must be finite."
                );
            }
            const double weight = area_scale * quadrature.weight;
            if (problem.compute_energy_balance) {
                neumann_heat_rate += static_cast<long double>(flux) * weight;
            }
            for (std::size_t i = 0; i < 3U; ++i) {
                const std::size_t global_i = face[i];
                const double boundary_load =
                    flux * quadrature.shape[i] * weight;
                rhs[global_i] -= boundary_load;
                nullspace.record_load(global_i, boundary_load);
                if (!std::isfinite(rhs[global_i])) {
                    throw std::overflow_error(
                        "3D Neumann load-vector accumulation overflowed."
                    );
                }
            }
        }
    };

    const auto integrate_robin_face = [&](
        const BoundaryFace3D& face,
        const ScalarField3D& transfer_coefficient,
        const ScalarField3D& ambient_value
    ) {
        const double area_scale = detail::face_jacobian(problem.mesh, face);
        for (const auto& quadrature : triangle_face_degree_two_quadrature) {
            const auto point = detail::simplex_point(problem.mesh, face, quadrature.shape);
            const double transfer = transfer_coefficient(
                point.x, point.y, point.z
            );
            const double ambient = ambient_value(
                point.x, point.y, point.z
            );
            if (!std::isfinite(transfer) || transfer < 0.0) {
                throw std::invalid_argument(
                    "3D Robin transfer coefficient must be finite and nonnegative."
                );
            }
            if (transfer > 0.0) {
                nullspace.anchor(face[0]);
            }
            if (!std::isfinite(ambient)) {
                throw std::invalid_argument(
                    "3D Robin ambient value must be finite."
                );
            }
            const double ambient_load = transfer * ambient;
            if (!std::isfinite(ambient_load)) {
                throw std::overflow_error(
                    "3D Robin ambient load is non-finite."
                );
            }

            const double weight = area_scale * quadrature.weight;
            if (problem.compute_energy_balance) {
                robin_energy_samples.push_back({
                    face, quadrature.shape, transfer, ambient, weight
                });
            }
            for (std::size_t i = 0; i < 3U; ++i) {
                const std::size_t global_i = face[i];
                const double boundary_load =
                    ambient_load * quadrature.shape[i] * weight;
                rhs[global_i] += boundary_load;
                nullspace.record_load(global_i, boundary_load);
                if (!std::isfinite(rhs[global_i])) {
                    throw std::overflow_error(
                        "3D Robin load-vector accumulation overflowed."
                    );
                }
                for (std::size_t j = 0; j < 3U; ++j) {
                    const double value =
                        transfer * quadrature.shape[i] *
                        quadrature.shape[j] * weight;
                    if (!std::isfinite(value)) {
                        throw std::overflow_error(
                            "3D Robin matrix entry is non-finite."
                        );
                    }
                    entries.push_back({global_i, face[j], value});
                }
            }
        }
    };

    for (const NeumannBoundaryFace3D& condition : problem.neumann_faces) {
        integrate_neumann_face(condition.face, condition.outward_flux);
    }

    for (const RobinBoundaryFace3D& condition : problem.robin_faces) {
        integrate_robin_face(
            condition.face,
            condition.transfer_coefficient,
            condition.ambient_value
        );
    }

    // Capture only true physical Dirichlet rows before temporary nullspace
    // gauges and symmetric elimination modify the assembled linear system.
    // Retaining just constrained rows avoids copying the full triplet matrix.
    std::vector<linalg::Triplet> dirichlet_rows;
    long double dirichlet_rhs_sum = 0.0L;
    if (problem.compute_energy_balance) {
        for (std::size_t node = 0; node < node_count; ++node) {
            if (prescribed[node]) dirichlet_rhs_sum += rhs[node];
        }
        for (const auto& entry : entries) {
            if (prescribed[entry.row]) dirichlet_rows.push_back(entry);
        }
    }

    nullspace.impose_gauges(rhs, prescribed);

    auto solver = detail::solve_constrained_diffusion(
        std::move(entries), std::move(rhs), prescribed, solver_options);
    auto values = nullspace.physical_values(solver.solution);
    std::optional<EnergyBalance3D> balance;
    if (problem.compute_energy_balance) {
        // A*u - rhs at a Dirichlet node is minus its physical
        // outward heat flux: integration by parts gives the sign.
        long double dirichlet_heat_rate = dirichlet_rhs_sum;
        for (const auto& entry : dirichlet_rows) {
            dirichlet_heat_rate -= static_cast<long double>(entry.value) *
                values[entry.column];
        }
        long double robin_heat_rate = 0.0L;
        for (const auto& sample : robin_energy_samples) {
            long double temperature = 0.0L;
            for (std::size_t i = 0; i < 3U; ++i) {
                temperature += static_cast<long double>(sample.shape[i]) *
                    values[sample.face[i]];
            }
            robin_heat_rate += static_cast<long double>(sample.transfer) *
                (temperature - sample.ambient) * sample.weight;
        }
        const long double net = dirichlet_heat_rate +
            neumann_heat_rate + robin_heat_rate;
        const long double error = net - generated_heat_rate;
        const long double scale = std::max({
            std::abs(generated_heat_rate), std::abs(dirichlet_heat_rate),
            std::abs(neumann_heat_rate), std::abs(robin_heat_rate)
        });
        balance = EnergyBalance3D{
            finite_energy_rate(generated_heat_rate),
            finite_energy_rate(dirichlet_heat_rate),
            finite_energy_rate(neumann_heat_rate),
            finite_energy_rate(robin_heat_rate),
            finite_energy_rate(net),
            finite_energy_rate(error),
            finite_energy_rate(scale == 0.0L ? 0.0L : std::abs(error) / scale)
        };
    }
    return {std::move(values), std::move(solver), std::move(balance)};
}

} // namespace nova::fem
