#include "nova/fem/mesh_transfer_2d.hpp"

#include <algorithm>
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

} // namespace

RefinedTransientHeat2DProblem
refine_transient_heat_problem_adaptively_2d(
    const TransientHeat2DProblem& problem,
    const H1TriangleOrder order,
    const std::span<const std::size_t> marked_elements
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
    result.refinement = refine_triangles_adaptively_2d(
        problem.mesh, marked_elements, assigned_edges
    );
    result.problem = problem;
    result.problem.mesh = result.refinement.mesh;

    if (!problem.initial_dof_values.empty()) {
        result.problem.initial_dof_values = prolongate_h1_refined_2d(
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

    const auto remap_edge = [&](const BoundaryEdge2D& edge) {
        std::vector<BoundaryEdge2D> segments;
        const auto midpoint = midpoint_ids.find(
            canonical_edge(edge[0], edge[1])
        );
        if (midpoint == midpoint_ids.end()) {
            segments.push_back(edge);
        } else {
            segments.push_back({edge[0], midpoint->second});
            segments.push_back({midpoint->second, edge[1]});
        }
        return segments;
    };

    result.problem.dirichlet_edges.clear();
    for (const auto& condition : problem.dirichlet_edges) {
        for (const BoundaryEdge2D& edge : remap_edge(condition.edge)) {
            auto transferred = condition;
            transferred.edge = edge;
            result.problem.dirichlet_edges.push_back(std::move(transferred));
        }
    }

    result.problem.neumann_edges.clear();
    for (const auto& condition : problem.neumann_edges) {
        for (const BoundaryEdge2D& edge : remap_edge(condition.edge)) {
            auto transferred = condition;
            transferred.edge = edge;
            result.problem.neumann_edges.push_back(std::move(transferred));
        }
    }

    result.problem.robin_edges.clear();
    for (const auto& condition : problem.robin_edges) {
        for (const BoundaryEdge2D& edge : remap_edge(condition.edge)) {
            auto transferred = condition;
            transferred.edge = edge;
            result.problem.robin_edges.push_back(std::move(transferred));
        }
    }

    return result;
}

} // namespace nova::fem
