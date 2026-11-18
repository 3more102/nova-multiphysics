#include "nova/fem/mesh_refinement_2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {

using EdgeKey = std::pair<std::size_t, std::size_t>;

[[nodiscard]] EdgeKey canonical_edge(
    const std::size_t first,
    const std::size_t second
) {
    return std::minmax(first, second);
}

[[nodiscard]] std::map<EdgeKey, std::size_t> edge_incidence(
    const TriangleMesh2D& mesh
) {
    std::map<EdgeKey, std::size_t> incidence;
    for (const Triangle2D& triangle : mesh.elements) {
        ++incidence[canonical_edge(triangle[0], triangle[1])];
        ++incidence[canonical_edge(triangle[1], triangle[2])];
        ++incidence[canonical_edge(triangle[2], triangle[0])];
    }
    return incidence;
}

} // namespace

RefinedTriangleMesh2D refine_triangles_adaptively_2d(
    const TriangleMesh2D& mesh,
    const std::span<const std::size_t> marked_elements,
    const std::vector<BoundaryEdge2D>& boundary_edges
) {
    mesh.validate();
    const auto incidence = edge_incidence(mesh);

    std::set<EdgeKey> unique_boundaries;
    for (const BoundaryEdge2D& edge : boundary_edges) {
        if (edge[0] >= mesh.nodes.size() || edge[1] >= mesh.nodes.size()) {
            throw std::out_of_range("Adaptive-refinement boundary references a missing node.");
        }
        if (edge[0] == edge[1]) {
            throw std::invalid_argument("Adaptive-refinement boundary repeats a node.");
        }
        const EdgeKey key = canonical_edge(edge[0], edge[1]);
        const auto found = incidence.find(key);
        if (found == incidence.end() || found->second != 1U) {
            throw std::invalid_argument(
                "Adaptive-refinement boundary edge must be exterior."
            );
        }
        if (!unique_boundaries.insert(key).second) {
            throw std::invalid_argument(
                "Adaptive-refinement boundary edge is duplicated."
            );
        }
    }

    std::vector<bool> marked(mesh.elements.size(), false);
    std::set<EdgeKey> split_edges;
    for (const std::size_t element : marked_elements) {
        if (element >= mesh.elements.size()) {
            throw std::out_of_range("Adaptive refinement references a missing element.");
        }
        if (marked[element]) {
            throw std::invalid_argument("Adaptive refinement contains a duplicate element mark.");
        }
        marked[element] = true;
        const Triangle2D& triangle = mesh.elements[element];
        split_edges.insert(canonical_edge(triangle[0], triangle[1]));
        split_edges.insert(canonical_edge(triangle[1], triangle[2]));
        split_edges.insert(canonical_edge(triangle[2], triangle[0]));
    }

    RefinedTriangleMesh2D result;
    result.coarse_node_count = mesh.nodes.size();
    result.mesh.nodes = mesh.nodes;

    if (marked_elements.empty()) {
        result.mesh.elements = mesh.elements;
        result.parent_elements.resize(mesh.elements.size());
        std::iota(
            result.parent_elements.begin(),
            result.parent_elements.end(),
            std::size_t{0}
        );
        result.boundary_edges = boundary_edges;
        return result;
    }

    if (split_edges.size() > result.mesh.nodes.max_size() - mesh.nodes.size() ||
        mesh.elements.size() > result.mesh.elements.max_size() / 4U ||
        mesh.elements.size() > result.parent_elements.max_size() / 4U ||
        boundary_edges.size() > result.boundary_edges.max_size() / 2U) {
        throw std::length_error("Adaptive triangle refinement exceeds container capacity.");
    }

    result.mesh.nodes.reserve(mesh.nodes.size() + split_edges.size());
    result.midpoint_parent_edges.reserve(split_edges.size());
    result.mesh.elements.reserve(4U * mesh.elements.size());
    result.parent_elements.reserve(4U * mesh.elements.size());
    result.boundary_edges.reserve(2U * boundary_edges.size());

    std::map<EdgeKey, std::size_t> midpoint_ids;
    for (const EdgeKey& edge : split_edges) {
        const Point2D& first = mesh.nodes[edge.first];
        const Point2D& second = mesh.nodes[edge.second];
        const double x = 0.5 * first.x + 0.5 * second.x;
        const double y = 0.5 * first.y + 0.5 * second.y;
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw std::overflow_error("Adaptive-refinement midpoint is non-finite.");
        }
        const std::size_t midpoint = result.mesh.nodes.size();
        midpoint_ids.emplace(edge, midpoint);
        result.mesh.nodes.push_back({x, y});
        result.midpoint_parent_edges.push_back({edge.first, edge.second});
    }

    const auto append_child = [&](const std::size_t parent, const Triangle2D& child) {
        result.mesh.elements.push_back(child);
        result.parent_elements.push_back(parent);
    };

    for (std::size_t parent = 0; parent < mesh.elements.size(); ++parent) {
        const Triangle2D& triangle = mesh.elements[parent];
        const std::size_t a = triangle[0];
        const std::size_t b = triangle[1];
        const std::size_t c = triangle[2];
        const EdgeKey edge_ab = canonical_edge(a, b);
        const EdgeKey edge_bc = canonical_edge(b, c);
        const EdgeKey edge_ca = canonical_edge(c, a);
        const bool split_ab = split_edges.contains(edge_ab);
        const bool split_bc = split_edges.contains(edge_bc);
        const bool split_ca = split_edges.contains(edge_ca);
        const unsigned int split_count =
            static_cast<unsigned int>(split_ab) +
            static_cast<unsigned int>(split_bc) +
            static_cast<unsigned int>(split_ca);

        if (split_count == 0U) {
            append_child(parent, triangle);
            continue;
        }

        const std::size_t mab = split_ab ? midpoint_ids.at(edge_ab) : 0U;
        const std::size_t mbc = split_bc ? midpoint_ids.at(edge_bc) : 0U;
        const std::size_t mca = split_ca ? midpoint_ids.at(edge_ca) : 0U;

        if (split_count == 1U) {
            if (split_ab) {
                append_child(parent, {a, mab, c});
                append_child(parent, {mab, b, c});
            } else if (split_bc) {
                append_child(parent, {b, mbc, a});
                append_child(parent, {mbc, c, a});
            } else {
                append_child(parent, {c, mca, b});
                append_child(parent, {mca, a, b});
            }
            continue;
        }

        if (split_count == 2U) {
            if (split_ab && split_bc) {
                append_child(parent, {b, mbc, mab});
                append_child(parent, {a, mab, c});
                append_child(parent, {mab, mbc, c});
            } else if (split_bc && split_ca) {
                append_child(parent, {c, mca, mbc});
                append_child(parent, {b, mbc, a});
                append_child(parent, {mbc, mca, a});
            } else {
                append_child(parent, {a, mab, mca});
                append_child(parent, {c, mca, b});
                append_child(parent, {mca, mab, b});
            }
            continue;
        }

        append_child(parent, {a, mab, mca});
        append_child(parent, {mab, b, mbc});
        append_child(parent, {mca, mbc, c});
        append_child(parent, {mab, mbc, mca});
    }

    // Keep each caller-supplied boundary edge in its original direction.
    // It appears once when untouched and as two segments when bisected.
    for (const BoundaryEdge2D& edge : boundary_edges) {
        const auto midpoint = midpoint_ids.find(
            canonical_edge(edge[0], edge[1])
        );
        if (midpoint == midpoint_ids.end()) {
            result.boundary_edges.push_back(edge);
        } else {
            result.boundary_edges.push_back({edge[0], midpoint->second});
            result.boundary_edges.push_back({midpoint->second, edge[1]});
        }
    }

    result.mesh.validate();
    return result;
}

} // namespace nova::fem
