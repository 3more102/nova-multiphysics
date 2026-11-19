#include "nova/fem/mesh_refinement_2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using EdgeKey = std::pair<std::size_t, std::size_t>;

[[nodiscard]] EdgeKey canonical_edge(std::size_t first, std::size_t second) {
    return std::minmax(first, second);
}

[[nodiscard]] double angle_between(
    const Point2D& origin, const Point2D& first, const Point2D& second
) {
    const double ax = first.x - origin.x;
    const double ay = first.y - origin.y;
    const double bx = second.x - origin.x;
    const double by = second.y - origin.y;
    return std::atan2(
        std::abs(ax * by - ay * bx),
        ax * bx + ay * by
    );
}

[[nodiscard]] std::map<EdgeKey, std::size_t> edge_incidence(
    const TriangleMesh2D& mesh
) {
    std::map<EdgeKey, std::size_t> counts;
    for (const Triangle2D& triangle : mesh.elements) {
        ++counts[canonical_edge(triangle[0], triangle[1])];
        ++counts[canonical_edge(triangle[1], triangle[2])];
        ++counts[canonical_edge(triangle[2], triangle[0])];
    }
    return counts;
}

} // namespace

TriangleMeshQuality2D evaluate_mesh_quality_2d(const TriangleMesh2D& mesh) {
    mesh.validate();
    TriangleMeshQuality2D result;
    result.elements.reserve(mesh.elements.size());
    result.minimum_angle_radians = std::numeric_limits<double>::infinity();
    result.minimum_mean_ratio_quality = 1.0;

    constexpr double sqrt_three = 1.7320508075688772935;
    for (const Triangle2D& triangle : mesh.elements) {
        const Point2D& a = mesh.nodes[triangle[0]];
        const Point2D& b = mesh.nodes[triangle[1]];
        const Point2D& c = mesh.nodes[triangle[2]];
        const double ab = std::hypot(b.x - a.x, b.y - a.y);
        const double bc = std::hypot(c.x - b.x, c.y - b.y);
        const double ca = std::hypot(a.x - c.x, a.y - c.y);
        const double det =
            (b.x - a.x) * (c.y - a.y) -
            (b.y - a.y) * (c.x - a.x);
        const double sum_squared = ab * ab + bc * bc + ca * ca;
        if (!std::isfinite(sum_squared) || sum_squared <= 0.0) {
            throw std::overflow_error("Triangle quality squared lengths overflowed.");
        }
        const double minimum_angle = std::min({
            angle_between(a, b, c),
            angle_between(b, c, a),
            angle_between(c, a, b)
        });
        const double quality = std::clamp(
            (2.0 * sqrt_three) * (std::abs(det) / sum_squared),
            0.0, 1.0
        );
        const double area = 0.5 * std::abs(det);
        result.elements.push_back({
            0.5 * det, minimum_angle, quality,
            std::max({ab, bc, ca})
        });
        result.minimum_angle_radians = std::min(
            result.minimum_angle_radians, minimum_angle
        );
        result.minimum_mean_ratio_quality = std::min(
            result.minimum_mean_ratio_quality, quality
        );
        result.maximum_edge_length = std::max(
            result.maximum_edge_length, std::max({ab, bc, ca})
        );
        result.total_area += area;
        if (!std::isfinite(result.total_area)) {
            throw std::overflow_error("Total triangle mesh area overflowed.");
        }
    }
    return result;
}

RefinedTriangleMesh2D refine_triangles_uniformly_2d(
    const TriangleMesh2D& mesh,
    const std::vector<BoundaryEdge2D>& boundary_edges
) {
    mesh.validate();
    const auto incidence = edge_incidence(mesh);

    std::set<EdgeKey> used_boundaries;
    for (const BoundaryEdge2D& edge : boundary_edges) {
        if (edge[0] >= mesh.nodes.size() || edge[1] >= mesh.nodes.size()) {
            throw std::out_of_range("Refinement boundary references a missing node.");
        }
        if (edge[0] == edge[1]) {
            throw std::invalid_argument("Refinement boundary edge repeats a node.");
        }
        const EdgeKey key = canonical_edge(edge[0], edge[1]);
        const auto it = incidence.find(key);
        if (it == incidence.end() || it->second != 1U) {
            throw std::invalid_argument(
                "Refinement boundary edge must be an exterior mesh edge."
            );
        }
        if (!used_boundaries.insert(key).second) {
            throw std::invalid_argument("Refinement boundary edge is duplicated.");
        }
    }

    RefinedTriangleMesh2D result;
    result.coarse_node_count = mesh.nodes.size();

    if (incidence.size() > result.mesh.nodes.max_size() - mesh.nodes.size() ||
        mesh.elements.size() > result.mesh.elements.max_size() / 4U ||
        mesh.elements.size() > result.parent_elements.max_size() / 4U ||
        boundary_edges.size() > result.boundary_edges.max_size() / 2U) {
        throw std::length_error("Uniform triangle refinement exceeds capacity.");
    }

    result.mesh.nodes = mesh.nodes;
    result.mesh.nodes.reserve(mesh.nodes.size() + incidence.size());
    result.midpoint_parent_edges.reserve(incidence.size());
    result.mesh.elements.reserve(4U * mesh.elements.size());
    result.parent_elements.reserve(4U * mesh.elements.size());
    result.boundary_edges.reserve(2U * boundary_edges.size());

    // Deterministic canonical edge traversal, independent of cell order.
    std::map<EdgeKey, std::size_t> midpoint_ids;
    for (const auto& [edge, count] : incidence) {
        static_cast<void>(count);
        const Point2D& a = mesh.nodes[edge.first];
        const Point2D& b = mesh.nodes[edge.second];
        const double x = 0.5 * a.x + 0.5 * b.x;
        const double y = 0.5 * a.y + 0.5 * b.y;
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw std::overflow_error("Refined mesh midpoint is non-finite.");
        }
        const std::size_t midpoint = result.mesh.nodes.size();
        midpoint_ids.emplace(edge, midpoint);
        result.mesh.nodes.push_back({x, y});
        result.midpoint_parent_edges.push_back({
            edge.first, edge.second
        });
    }

    for (std::size_t parent = 0; parent < mesh.elements.size(); ++parent) {
        const auto& t = mesh.elements[parent];
        const std::size_t a = t[0];
        const std::size_t b = t[1];
        const std::size_t c = t[2];
        const std::size_t ab = midpoint_ids.at(canonical_edge(a, b));
        const std::size_t bc = midpoint_ids.at(canonical_edge(b, c));
        const std::size_t ca = midpoint_ids.at(canonical_edge(c, a));
        result.mesh.elements.push_back({a, ab, ca});
        result.mesh.elements.push_back({ab, b, bc});
        result.mesh.elements.push_back({ca, bc, c});
        result.mesh.elements.push_back({ab, bc, ca});
        for (std::size_t child = 0; child < 4U; ++child) {
            result.parent_elements.push_back(parent);
        }
    }

    for (const BoundaryEdge2D& edge : boundary_edges) {
        const std::size_t midpoint =
            midpoint_ids.at(canonical_edge(edge[0], edge[1]));
        result.boundary_edges.push_back({edge[0], midpoint});
        result.boundary_edges.push_back({midpoint, edge[1]});
    }

    result.mesh.validate();
    return result;
}

std::vector<double> prolongate_p1_uniform_2d(
    const RefinedTriangleMesh2D& refined,
    const std::span<const double> coarse_nodal_values
) {
    if (coarse_nodal_values.size() != refined.coarse_node_count ||
        refined.mesh.nodes.size() !=
            refined.coarse_node_count + refined.midpoint_parent_edges.size()) {
        throw std::invalid_argument("Uniform P1 prolongation has inconsistent DOF sizes.");
    }
    std::vector<double> values(
        coarse_nodal_values.begin(), coarse_nodal_values.end()
    );
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Coarse P1 field must contain finite values.");
        }
    }
    values.reserve(refined.mesh.nodes.size());
    for (const BoundaryEdge2D& edge : refined.midpoint_parent_edges) {
        if (edge[0] >= refined.coarse_node_count ||
            edge[1] >= refined.coarse_node_count) {
            throw std::out_of_range("Prolongation references a missing coarse DOF.");
        }
        const double midpoint_value =
            0.5 * values[edge[0]] + 0.5 * values[edge[1]];
        if (!std::isfinite(midpoint_value)) {
            throw std::overflow_error("Prolongated P1 value is non-finite.");
        }
        values.push_back(midpoint_value);
    }
    return values;
}

} // namespace nova::fem
