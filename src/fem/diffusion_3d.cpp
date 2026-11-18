#include "nova/fem/diffusion_3d.hpp"

#include "p1_integration_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

class DisjointSet {
public:
    explicit DisjointSet(const std::size_t size)
        : parent_(size), rank_(size, 0U) {
        for (std::size_t node = 0; node < size; ++node) {
            parent_[node] = node;
        }
    }

    [[nodiscard]] std::size_t find(const std::size_t node) {
        if (parent_[node] != node) {
            parent_[node] = find(parent_[node]);
        }
        return parent_[node];
    }

    void unite(const std::size_t first, const std::size_t second) {
        std::size_t first_root = find(first);
        std::size_t second_root = find(second);
        if (first_root == second_root) {
            return;
        }
        if (rank_[first_root] < rank_[second_root]) {
            std::swap(first_root, second_root);
        }
        parent_[second_root] = first_root;
        if (rank_[first_root] == rank_[second_root]) {
            ++rank_[first_root];
        }
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<std::size_t> rank_;
};

struct MeshComponents {
    std::vector<std::size_t> node_component;
    std::vector<std::vector<std::size_t>> nodes;
};

[[nodiscard]] MeshComponents find_mesh_components(
    const TetrahedralMesh3D& mesh
) {
    DisjointSet sets(mesh.nodes.size());
    for (const Tetrahedron3D& cell : mesh.elements) {
        for (const std::size_t node : cell) {
            sets.unite(cell[0], node);
        }
    }

    MeshComponents result;
    result.node_component.resize(mesh.nodes.size());
    std::map<std::size_t, std::size_t> component_for_root;
    for (std::size_t node = 0; node < mesh.nodes.size(); ++node) {
        const std::size_t root = sets.find(node);
        const auto [entry, inserted] = component_for_root.emplace(
            root, component_for_root.size()
        );
        if (inserted) {
            result.nodes.emplace_back();
        }
        const std::size_t component = entry->second;
        result.node_component[node] = component;
        result.nodes[component].push_back(node);
    }
    return result;
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
    const MeshComponents components = find_mesh_components(problem.mesh);
    const std::size_t component_count = components.nodes.size();
    std::vector<bool> anchored_component(component_count, false);
    std::vector<long double> absolute_load_scale(component_count, 0.0L);
    std::vector<long double> nodal_volume_weights(node_count, 0.0L);

    if (!problem.diffusivity || !problem.source) {
        throw std::invalid_argument(
            "3D diffusivity and source functions are required."
        );
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

    for (std::size_t node = 0; node < node_count; ++node) {
        if (prescribed[node].has_value()) {
            anchored_component[components.node_component[node]] = true;
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
            nodal_volume_weights[node] += nodal_volume_weight;
        }

        std::array<std::array<double, 4>, 4> local_matrix{};
        std::array<double, 4> local_rhs{};

        for (const auto& quadrature : tetrahedron_degree_two_quadrature) {
            const auto point = detail::simplex_point(problem.mesh, cell, quadrature.shape);
            const double diffusivity = problem.diffusivity(point.x, point.y, point.z);
            const double source = problem.source(point.x, point.y, point.z);
            if (!std::isfinite(diffusivity) || diffusivity <= 0.0) {
                throw std::invalid_argument(
                    "3D diffusivity must be finite and strictly positive."
                );
            }
            if (!std::isfinite(source)) {
                throw std::invalid_argument("3D source must be finite.");
            }

            const double weight = absolute_determinant * quadrature.weight;
            const std::size_t component =
                components.node_component[cell[0]];
            for (std::size_t i = 0; i < local_dof_count; ++i) {
                const double source_load =
                    source * quadrature.shape[i] * weight;
                local_rhs[i] += source_load;
                absolute_load_scale[component] += std::abs(
                    static_cast<long double>(source_load)
                );
                for (std::size_t j = 0; j < local_dof_count; ++j) {
                    local_matrix[i][j] +=
                        diffusivity * detail::gradient_product(physical, i, j) * weight;
                }
            }
        }

        for (std::size_t i = 0; i < local_dof_count; ++i) {
            const std::size_t global_i = cell[i];
            rhs[global_i] += local_rhs[i];
            if (!std::isfinite(rhs[global_i])) {
                throw std::overflow_error(
                    "3D element load-vector accumulation overflowed."
                );
            }
            for (std::size_t j = 0; j < local_dof_count; ++j) {
                if (!std::isfinite(local_matrix[i][j])) {
                    throw std::overflow_error(
                        "3D element stiffness entry is non-finite."
                    );
                }
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
        const std::size_t component =
            components.node_component[face[0]];
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
            for (std::size_t i = 0; i < 3U; ++i) {
                const std::size_t global_i = face[i];
                const double boundary_load =
                    flux * quadrature.shape[i] * weight;
                rhs[global_i] -= boundary_load;
                absolute_load_scale[component] += std::abs(
                    static_cast<long double>(boundary_load)
                );
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
        const std::size_t component =
            components.node_component[face[0]];
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
                anchored_component[component] = true;
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
            for (std::size_t i = 0; i < 3U; ++i) {
                const std::size_t global_i = face[i];
                const double boundary_load =
                    ambient_load * quadrature.shape[i] * weight;
                rhs[global_i] += boundary_load;
                absolute_load_scale[component] += std::abs(
                    static_cast<long double>(boundary_load)
                );
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

    std::vector<std::size_t> pure_neumann_components;
    for (std::size_t component = 0;
         component < component_count;
         ++component) {
        if (anchored_component[component]) {
            continue;
        }

        long double net_load = 0.0L;
        for (const std::size_t node : components.nodes[component]) {
            net_load += static_cast<long double>(rhs[node]);
        }
        constexpr long double compatibility_relative_tolerance = 1.0e-10L;
        const long double compatibility_tolerance =
            compatibility_relative_tolerance * absolute_load_scale[component];
        if (std::abs(net_load) > compatibility_tolerance) {
            throw std::invalid_argument(
                "Pure-Neumann 3D component " + std::to_string(component) +
                " is incompatible: the integrated source must equal the "
                "integrated outward Neumann flux."
            );
        }

        pure_neumann_components.push_back(component);
        const std::size_t gauge_node = components.nodes[component].front();
        prescribed[gauge_node] = 0.0;
    }

    // Symmetric strong Dirichlet elimination preserves CG's SPD requirement.
    std::vector<linalg::Triplet> constrained_entries;
    if (node_count > constrained_entries.max_size() - entries.size()) {
        throw std::length_error("Constrained 3D finite-element matrix is too large.");
    }
    constrained_entries.reserve(entries.size() + node_count);
    for (const linalg::Triplet& entry : entries) {
        if (!prescribed[entry.row].has_value()) {
            if (prescribed[entry.column].has_value()) {
                rhs[entry.row] -=
                    entry.value * *prescribed[entry.column];
                if (!std::isfinite(rhs[entry.row])) {
                    throw std::overflow_error(
                        "3D Dirichlet right-hand-side correction overflowed."
                    );
                }
            } else {
                constrained_entries.push_back(entry);
            }
        }
    }
    for (std::size_t node = 0; node < node_count; ++node) {
        if (prescribed[node].has_value()) {
            rhs[node] = *prescribed[node];
            constrained_entries.push_back({node, node, 1.0});
        }
    }

    const linalg::CsrMatrix matrix = linalg::CsrMatrix::from_triplets(
        node_count, node_count, constrained_entries
    );
    linalg::SolverResult solver = linalg::solve_conjugate_gradient(
        matrix, std::span<const double>(rhs), solver_options
    );

    auto values = solver.solution;
    for (const std::size_t component : pure_neumann_components) {
        long double weighted_value_sum = 0.0L;
        long double component_volume = 0.0L;
        for (const std::size_t node : components.nodes[component]) {
            const long double weight = nodal_volume_weights[node];
            weighted_value_sum +=
                weight * static_cast<long double>(values[node]);
            component_volume += weight;
        }
        if (!std::isfinite(component_volume) || component_volume <= 0.0L) {
            throw std::logic_error(
                "Pure-Neumann 3D component has no finite positive volume."
            );
        }
        const long double weighted_mean =
            weighted_value_sum / component_volume;
        const double mean = static_cast<double>(weighted_mean);
        if (!std::isfinite(weighted_mean) || !std::isfinite(mean)) {
            throw std::overflow_error(
                "Pure-Neumann 3D component mean is non-finite."
            );
        }
        for (const std::size_t node : components.nodes[component]) {
            values[node] -= mean;
            if (!std::isfinite(values[node])) {
                throw std::overflow_error(
                    "Pure-Neumann 3D zero-mean normalization overflowed."
                );
            }
        }
    }

    return {std::move(values), std::move(solver)};
}

} // namespace nova::fem
