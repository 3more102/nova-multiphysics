#include "nova/fem/mesh_transfer_2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using EdgeKey = EdgeKey2D;

[[nodiscard]] EdgeKey canonical_edge(
    const std::size_t first,
    const std::size_t second
) {
    return std::minmax(first, second);
}

[[nodiscard]] std::array<double, 2> reference_coordinates_for_node(
    const TriangleMesh2D& coarse_mesh,
    const RefinedTriangleMesh2D& refined,
    const std::size_t parent_element,
    const std::size_t fine_node
) {
    const Triangle2D& triangle = coarse_mesh.elements.at(parent_element);
    if (fine_node < refined.coarse_node_count) {
        for (std::size_t local = 0; local < 3U; ++local) {
            if (triangle[local] == fine_node) {
                if (local == 0U) {
                    return {0.0, 0.0};
                }
                if (local == 1U) {
                    return {1.0, 0.0};
                }
                return {0.0, 1.0};
            }
        }
        throw std::invalid_argument(
            "A refined edge endpoint is not a vertex of its recorded parent."
        );
    }

    const std::size_t midpoint_index = fine_node - refined.coarse_node_count;
    if (midpoint_index >= refined.midpoint_parent_edges.size()) {
        throw std::out_of_range(
            "A refined edge endpoint has no coarse-edge midpoint record."
        );
    }
    const BoundaryEdge2D& parent_edge =
        refined.midpoint_parent_edges[midpoint_index];
    const auto matches = [&](const std::size_t first, const std::size_t second) {
        return canonical_edge(first, second) ==
            canonical_edge(parent_edge[0], parent_edge[1]);
    };

    if (matches(triangle[0], triangle[1])) {
        return {0.5, 0.0};
    }
    if (matches(triangle[1], triangle[2])) {
        return {0.5, 0.5};
    }
    if (matches(triangle[2], triangle[0])) {
        return {0.0, 0.5};
    }
    throw std::invalid_argument(
        "A refined midpoint is not on an edge of its recorded parent."
    );
}

[[nodiscard]] double evaluate_coarse_field(
    const H1DofManager2D& coarse_dofs,
    const H1TriangleOrder order,
    const std::span<const double> coarse_values,
    const std::size_t parent_element,
    const double r,
    const double s
) {
    const TriangleH1BasisEvaluation basis =
        evaluate_triangle_h1_basis(order, r, s);
    const H1ElementDofMap2D& local_dofs =
        coarse_dofs.element_dofs(parent_element);

    double value = 0.0;
    for (std::size_t local = 0; local < basis.count; ++local) {
        value += basis.values[local] * coarse_values[local_dofs[local]];
    }
    if (!std::isfinite(value)) {
        throw std::overflow_error("H1 mesh-transfer evaluation is non-finite.");
    }
    return value;
}

} // namespace

std::vector<double> prolongate_h1_refined_2d(
    const TriangleMesh2D& coarse_mesh,
    const RefinedTriangleMesh2D& refined,
    const H1TriangleOrder order,
    const std::span<const double> coarse_dof_values
) {
    coarse_mesh.validate();
    refined.mesh.validate();

    if (refined.coarse_node_count != coarse_mesh.nodes.size() ||
        refined.coarse_node_count > refined.mesh.nodes.size() ||
        refined.mesh.nodes.size() - refined.coarse_node_count !=
            refined.midpoint_parent_edges.size() ||
        refined.parent_elements.size() != refined.mesh.elements.size()) {
        throw std::invalid_argument("Refined mesh transfer metadata is inconsistent.");
    }

    for (const std::size_t parent : refined.parent_elements) {
        if (parent >= coarse_mesh.elements.size()) {
            throw std::invalid_argument(
                "Refinement parent-element map references a missing coarse cell."
            );
        }
    }

    for (std::size_t node = 0; node < coarse_mesh.nodes.size(); ++node) {
        const Point2D& coarse = coarse_mesh.nodes[node];
        const Point2D& fine = refined.mesh.nodes[node];
        if (coarse.x != fine.x || coarse.y != fine.y) {
            throw std::invalid_argument(
                "Refined mesh does not preserve the original coarse vertices."
            );
        }
    }

    const H1FunctionSpace2D coarse_space(coarse_mesh, order);
    const H1FunctionSpace2D fine_space(refined.mesh, order);
    const H1DofManager2D& coarse_dofs = coarse_space.dofs();
    const H1DofManager2D& fine_dofs = fine_space.dofs();
    if (coarse_dof_values.size() != coarse_space.dimension()) {
        throw std::invalid_argument("Coarse H1 state has the wrong DOF count.");
    }
    for (const double value : coarse_dof_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Coarse H1 state must contain finite values.");
        }
    }

    std::vector<double> fine_values(fine_space.dimension(), 0.0);
    std::copy(
        coarse_dof_values.begin(),
        coarse_dof_values.begin() +
            static_cast<std::ptrdiff_t>(refined.coarse_node_count),
        fine_values.begin()
    );

    // Each newly created mesh vertex is the midpoint of one unique coarse
    // edge. P1 uses linear interpolation; P2 reuses that edge's exact nodal DOF.
    for (std::size_t edge_index = 0;
         edge_index < refined.midpoint_parent_edges.size();
         ++edge_index) {
        const BoundaryEdge2D& edge =
            refined.midpoint_parent_edges[edge_index];
        if (edge[0] >= refined.coarse_node_count ||
            edge[1] >= refined.coarse_node_count ||
            edge[0] == edge[1]) {
            throw std::invalid_argument(
                "Refinement midpoint references invalid coarse vertices."
            );
        }

        const std::size_t new_vertex = refined.coarse_node_count + edge_index;
        if (order == H1TriangleOrder::P1) {
            fine_values[new_vertex] =
                0.5 * coarse_dof_values[edge[0]] +
                0.5 * coarse_dof_values[edge[1]];
        } else {
            const auto coarse_edge_dof =
                coarse_dofs.edge_dof(edge[0], edge[1]);
            if (!coarse_edge_dof.has_value()) {
                throw std::invalid_argument(
                    "Refinement midpoint is not an edge of the coarse mesh."
                );
            }
            fine_values[new_vertex] = coarse_dof_values[*coarse_edge_dof];
        }
        if (!std::isfinite(fine_values[new_vertex])) {
            throw std::overflow_error("Transferred H1 vertex value is non-finite.");
        }
    }

    if (order == H1TriangleOrder::P1) {
        return fine_values;
    }

    // A fine P2 edge DOF is evaluated from its coarse parent polynomial.
    // Edges on coarse interfaces have two parents; continuity requires both
    // traces to agree, which is checked before accepting the transferred value.
    std::map<EdgeKey, std::vector<std::size_t>> edge_parents;
    for (std::size_t child = 0; child < refined.mesh.elements.size(); ++child) {
        const std::size_t parent = refined.parent_elements[child];
        if (parent >= coarse_mesh.elements.size()) {
            throw std::invalid_argument(
                "Refinement parent-element map references a missing coarse cell."
            );
        }
        const Triangle2D& triangle = refined.mesh.elements[child];
        const std::array<EdgeKey, 3> edges{{
            canonical_edge(triangle[0], triangle[1]),
            canonical_edge(triangle[1], triangle[2]),
            canonical_edge(triangle[2], triangle[0])
        }};
        for (const EdgeKey& edge : edges) {
            auto& parents = edge_parents[edge];
            if (std::find(parents.begin(), parents.end(), parent) == parents.end()) {
                parents.push_back(parent);
            }
        }
    }

    for (const EdgeKey& edge : fine_dofs.edge_keys()) {
        const auto parents = edge_parents.find(edge);
        if (parents == edge_parents.end() || parents->second.empty()) {
            throw std::invalid_argument(
                "A fine H1 edge has no parent-element transfer mapping."
            );
        }

        double transferred_value = 0.0;
        bool first_parent = true;
        for (const std::size_t parent : parents->second) {
            const auto first_reference = reference_coordinates_for_node(
                coarse_mesh, refined, parent, edge.first
            );
            const auto second_reference = reference_coordinates_for_node(
                coarse_mesh, refined, parent, edge.second
            );
            const double midpoint_r =
                0.5 * first_reference[0] + 0.5 * second_reference[0];
            const double midpoint_s =
                0.5 * first_reference[1] + 0.5 * second_reference[1];
            const double candidate = evaluate_coarse_field(
                coarse_dofs, order, coarse_dof_values,
                parent, midpoint_r, midpoint_s
            );
            if (first_parent) {
                transferred_value = candidate;
                first_parent = false;
            } else {
                const double scale = std::max({
                    1.0, std::abs(transferred_value), std::abs(candidate)
                });
                if (std::abs(candidate - transferred_value) >
                    1.0e-10 * scale) {
                    throw std::invalid_argument(
                        "Coarse H1 traces disagree across a refined shared edge."
                    );
                }
            }
        }

        const auto fine_edge_dof = fine_dofs.edge_dof(edge.first, edge.second);
        if (!fine_edge_dof.has_value()) {
            throw std::logic_error("Fine P2 edge DOF numbering is incomplete.");
        }
        fine_values[*fine_edge_dof] = transferred_value;
    }

    for (const double value : fine_values) {
        if (!std::isfinite(value)) {
            throw std::overflow_error("Transferred P2 state is non-finite.");
        }
    }
    return fine_values;
}


std::vector<double> prolongate_h1_uniform_2d(
    const TriangleMesh2D& coarse_mesh,
    const RefinedTriangleMesh2D& refined,
    const H1TriangleOrder order,
    const std::span<const double> coarse_dof_values
) {
    return prolongate_h1_refined_2d(
        coarse_mesh, refined, order, coarse_dof_values
    );
}

RefinedTransientHeat2DProblem
refine_transient_heat_problem_uniformly_2d(
    const TransientHeat2DProblem& problem,
    const H1TriangleOrder order
) {
    problem.mesh.validate();

    std::vector<BoundaryEdge2D> assigned_edges;
    std::set<EdgeKey> unique_edges;
    const auto collect_edge = [&](const BoundaryEdge2D& edge) {
        if (!unique_edges.insert(canonical_edge(edge[0], edge[1])).second) {
            throw std::invalid_argument(
                "A transient heat boundary edge is assigned more than once."
            );
        }
        assigned_edges.push_back(edge);
    };
    for (const auto& condition : problem.dirichlet_edges) {
        collect_edge(condition.edge);
    }
    for (const auto& condition : problem.neumann_edges) {
        collect_edge(condition.edge);
    }
    for (const auto& condition : problem.robin_edges) {
        collect_edge(condition.edge);
    }

    RefinedTransientHeat2DProblem result;
    result.refinement = refine_triangles_uniformly_2d(
        problem.mesh, assigned_edges
    );
    result.problem = problem;
    result.problem.mesh = result.refinement.mesh;

    if (!problem.initial_dof_values.empty()) {
        result.problem.initial_dof_values = prolongate_h1_uniform_2d(
            problem.mesh, result.refinement, order,
            problem.initial_dof_values
        );
    }

    std::map<EdgeKey, std::size_t> midpoint_ids;
    for (std::size_t i = 0;
         i < result.refinement.midpoint_parent_edges.size();
         ++i) {
        const BoundaryEdge2D& edge =
            result.refinement.midpoint_parent_edges[i];
        midpoint_ids.emplace(
            canonical_edge(edge[0], edge[1]),
            result.refinement.coarse_node_count + i
        );
    }

    const auto split_edge = [&](const BoundaryEdge2D& edge) {
        const auto midpoint = midpoint_ids.find(
            canonical_edge(edge[0], edge[1])
        );
        if (midpoint == midpoint_ids.end()) {
            throw std::logic_error(
                "A prescribed boundary edge has no refined midpoint."
            );
        }
        return std::array<BoundaryEdge2D, 2>{{
            BoundaryEdge2D{edge[0], midpoint->second},
            BoundaryEdge2D{midpoint->second, edge[1]}
        }};
    };

    result.problem.dirichlet_edges.clear();
    for (const auto& condition : problem.dirichlet_edges) {
        const auto segments = split_edge(condition.edge);
        auto first = condition;
        auto second = condition;
        first.edge = segments[0];
        second.edge = segments[1];
        result.problem.dirichlet_edges.push_back(std::move(first));
        result.problem.dirichlet_edges.push_back(std::move(second));
    }

    result.problem.neumann_edges.clear();
    for (const auto& condition : problem.neumann_edges) {
        const auto segments = split_edge(condition.edge);
        auto first = condition;
        auto second = condition;
        first.edge = segments[0];
        second.edge = segments[1];
        result.problem.neumann_edges.push_back(std::move(first));
        result.problem.neumann_edges.push_back(std::move(second));
    }

    result.problem.robin_edges.clear();
    for (const auto& condition : problem.robin_edges) {
        const auto segments = split_edge(condition.edge);
        auto first = condition;
        auto second = condition;
        first.edge = segments[0];
        second.edge = segments[1];
        result.problem.robin_edges.push_back(std::move(first));
        result.problem.robin_edges.push_back(std::move(second));
    }

    return result;
}

} // namespace nova::fem
