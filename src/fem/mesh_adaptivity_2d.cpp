#include "nova/fem/mesh_adaptivity_2d.hpp"

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

using EdgeKey = EdgeKey2D;

[[nodiscard]] EdgeKey canonical_edge(
    const std::size_t first, const std::size_t second
) {
    return std::minmax(first, second);
}

[[nodiscard]] double edge_length(
    const TriangleMesh2D& mesh, const EdgeKey& edge
) {
    const Point2D& a = mesh.nodes[edge.first];
    const Point2D& b = mesh.nodes[edge.second];
    return std::hypot(b.x - a.x, b.y - a.y);
}

struct ElementGradient2D {
    double x{};
    double y{};
    double area{};
    double diameter{};
};

[[nodiscard]] ElementGradient2D gradient_at_centroid(
    const TriangleMesh2D& mesh,
    const H1DofManager2D& dofs,
    const std::span<const double> values,
    const std::size_t element,
    const H1TriangleOrder order
) {
    const Triangle2D& triangle = mesh.elements[element];
    const Point2D& p0 = mesh.nodes[triangle[0]];
    const Point2D& p1 = mesh.nodes[triangle[1]];
    const Point2D& p2 = mesh.nodes[triangle[2]];

    const double j00 = p1.x - p0.x;
    const double j01 = p2.x - p0.x;
    const double j10 = p1.y - p0.y;
    const double j11 = p2.y - p0.y;
    const double determinant = j00 * j11 - j01 * j10;
    if (!std::isfinite(determinant) || determinant == 0.0) {
        throw std::invalid_argument("Cannot recover a gradient on a degenerate triangle.");
    }

    const TriangleH1BasisEvaluation basis =
        evaluate_triangle_h1_basis(order, 1.0 / 3.0, 1.0 / 3.0);
    const H1ElementDofMap2D& local = dofs.element_dofs(element);
    double grad_r = 0.0;
    double grad_s = 0.0;
    for (std::size_t i = 0; i < basis.count; ++i) {
        grad_r += values[local[i]] * basis.reference_gradients[i][0];
        grad_s += values[local[i]] * basis.reference_gradients[i][1];
    }

    const double grad_x = (j11 * grad_r - j10 * grad_s) / determinant;
    const double grad_y = (-j01 * grad_r + j00 * grad_s) / determinant;
    const double area = 0.5 * std::abs(determinant);
    const double diameter = std::max({
        edge_length(mesh, canonical_edge(triangle[0], triangle[1])),
        edge_length(mesh, canonical_edge(triangle[1], triangle[2])),
        edge_length(mesh, canonical_edge(triangle[2], triangle[0]))
    });

    if (!std::isfinite(grad_x) || !std::isfinite(grad_y) ||
        !std::isfinite(area) || !std::isfinite(diameter) ||
        area <= 0.0 || diameter <= 0.0) {
        throw std::overflow_error("Gradient recovery produced non-finite geometry or values.");
    }
    return {grad_x, grad_y, area, diameter};
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

std::vector<double> estimate_h1_gradient_recovery_indicators_2d(
    const TriangleMesh2D& mesh,
    const H1TriangleOrder order,
    const std::span<const double> dof_values
) {
    mesh.validate();
    const H1FunctionSpace2D space(mesh, order);
    if (dof_values.size() != space.dimension()) {
        throw std::invalid_argument("Gradient recovery received the wrong H1 DOF count.");
    }
    for (const double value : dof_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Gradient recovery requires finite H1 values.");
        }
    }

    const H1DofManager2D& dofs = space.dofs();
    std::vector<ElementGradient2D> gradients;
    gradients.reserve(mesh.elements.size());
    std::vector<double> recovered_x(mesh.nodes.size(), 0.0);
    std::vector<double> recovered_y(mesh.nodes.size(), 0.0);
    std::vector<double> weights(mesh.nodes.size(), 0.0);

    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const ElementGradient2D gradient = gradient_at_centroid(
            mesh, dofs, dof_values, element, order
        );
        gradients.push_back(gradient);
        for (const std::size_t vertex : mesh.elements[element]) {
            recovered_x[vertex] += gradient.area * gradient.x;
            recovered_y[vertex] += gradient.area * gradient.y;
            weights[vertex] += gradient.area;
        }
    }

    for (std::size_t vertex = 0; vertex < mesh.nodes.size(); ++vertex) {
        if (weights[vertex] > 0.0) {
            recovered_x[vertex] /= weights[vertex];
            recovered_y[vertex] /= weights[vertex];
        }
        if (!std::isfinite(recovered_x[vertex]) ||
            !std::isfinite(recovered_y[vertex])) {
            throw std::overflow_error("Recovered nodal gradient is non-finite.");
        }
    }

    std::vector<double> indicators(mesh.elements.size(), 0.0);
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const Triangle2D& triangle = mesh.elements[element];
        const double recovered_element_x = (
            recovered_x[triangle[0]] +
            recovered_x[triangle[1]] +
            recovered_x[triangle[2]]
        ) / 3.0;
        const double recovered_element_y = (
            recovered_y[triangle[0]] +
            recovered_y[triangle[1]] +
            recovered_y[triangle[2]]
        ) / 3.0;
        const double difference = std::hypot(
            recovered_element_x - gradients[element].x,
            recovered_element_y - gradients[element].y
        );
        const double indicator =
            gradients[element].diameter *
            std::sqrt(gradients[element].area) * difference;
        if (!std::isfinite(indicator)) {
            throw std::overflow_error("Gradient-recovery indicator is non-finite.");
        }
        indicators[element] = indicator;
    }
    return indicators;
}

std::vector<std::size_t> mark_elements_dorfler_2d(
    const std::span<const double> indicators,
    const double theta
) {
    if (!std::isfinite(theta) || theta <= 0.0 || theta > 1.0) {
        throw std::invalid_argument("Dorfler theta must lie in (0,1].");
    }

    std::vector<std::pair<long double, std::size_t>> ranked;
    ranked.reserve(indicators.size());
    long double total_squared = 0.0L;
    for (std::size_t element = 0; element < indicators.size(); ++element) {
        const double indicator = indicators[element];
        if (!std::isfinite(indicator) || indicator < 0.0) {
            throw std::invalid_argument(
                "Dorfler indicators must be finite and non-negative."
            );
        }
        const long double squared =
            static_cast<long double>(indicator) *
            static_cast<long double>(indicator);
        ranked.emplace_back(squared, element);
        total_squared += squared;
    }
    if (!std::isfinite(total_squared)) {
        throw std::overflow_error("Dorfler indicator norm overflowed.");
    }
    if (total_squared == 0.0L) {
        return {};
    }

    std::sort(
        ranked.begin(), ranked.end(),
        [](const auto& left, const auto& right) {
            if (left.first != right.first) {
                return left.first > right.first;
            }
            return left.second < right.second;
        }
    );

    const long double target =
        static_cast<long double>(theta) * total_squared;
    long double accumulated = 0.0L;
    std::vector<std::size_t> marked;
    for (const auto& [squared, element] : ranked) {
        if (squared == 0.0L) {
            break;
        }
        marked.push_back(element);
        accumulated += squared;
        if (accumulated >= target) {
            break;
        }
    }
    return marked;
}

RefinedTriangleMesh2D refine_triangles_by_marked_elements_2d(
    const TriangleMesh2D& mesh,
    const std::vector<BoundaryEdge2D>& boundary_edges,
    const std::span<const std::size_t> marked_elements
) {
    mesh.validate();
    const auto incidence = edge_incidence(mesh);

    std::set<std::size_t> unique_elements;
    std::set<EdgeKey> split_edges;
    for (const std::size_t element : marked_elements) {
        if (element >= mesh.elements.size()) {
            throw std::out_of_range("Adaptive refinement references a missing element.");
        }
        if (!unique_elements.insert(element).second) {
            throw std::invalid_argument("Adaptive refinement element is duplicated.");
        }
        const Triangle2D& triangle = mesh.elements[element];
        const std::array<EdgeKey, 3> edges{{
            canonical_edge(triangle[0], triangle[1]),
            canonical_edge(triangle[1], triangle[2]),
            canonical_edge(triangle[2], triangle[0])
        }};
        std::size_t longest = 0U;
        double longest_length = edge_length(mesh, edges[0]);
        for (std::size_t i = 1U; i < edges.size(); ++i) {
            const double length = edge_length(mesh, edges[i]);
            if (length > longest_length) {
                longest = i;
                longest_length = length;
            }
        }
        split_edges.insert(edges[longest]);
    }

    std::set<EdgeKey> unique_boundaries;
    for (const BoundaryEdge2D& edge : boundary_edges) {
        if (edge[0] >= mesh.nodes.size() || edge[1] >= mesh.nodes.size()) {
            throw std::out_of_range("Adaptive boundary references a missing node.");
        }
        if (edge[0] == edge[1]) {
            throw std::invalid_argument("Adaptive boundary edge repeats a node.");
        }
        const EdgeKey key = canonical_edge(edge[0], edge[1]);
        const auto found = incidence.find(key);
        if (found == incidence.end() || found->second != 1U) {
            throw std::invalid_argument(
                "Adaptive boundary edge must be an exterior mesh edge."
            );
        }
        if (!unique_boundaries.insert(key).second) {
            throw std::invalid_argument("Adaptive boundary edge is duplicated.");
        }
    }

    RefinedTriangleMesh2D result;
    result.coarse_node_count = mesh.nodes.size();
    if (split_edges.size() > result.mesh.nodes.max_size() - mesh.nodes.size() ||
        mesh.elements.size() > result.mesh.elements.max_size() / 4U ||
        mesh.elements.size() > result.parent_elements.max_size() / 4U) {
        throw std::length_error("Adaptive triangle refinement exceeds capacity.");
    }

    result.mesh.nodes = mesh.nodes;
    result.mesh.nodes.reserve(mesh.nodes.size() + split_edges.size());
    result.midpoint_parent_edges.reserve(split_edges.size());
    result.mesh.elements.reserve(4U * mesh.elements.size());
    result.parent_elements.reserve(4U * mesh.elements.size());

    std::map<EdgeKey, std::size_t> midpoint_ids;
    for (const EdgeKey& edge : split_edges) {
        const Point2D& a = mesh.nodes[edge.first];
        const Point2D& b = mesh.nodes[edge.second];
        const double x = 0.5 * a.x + 0.5 * b.x;
        const double y = 0.5 * a.y + 0.5 * b.y;
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw std::overflow_error("Adaptive midpoint is non-finite.");
        }
        const std::size_t midpoint = result.mesh.nodes.size();
        midpoint_ids.emplace(edge, midpoint);
        result.mesh.nodes.push_back({x, y});
        result.midpoint_parent_edges.push_back({edge.first, edge.second});
    }

    const auto midpoint = [&](const EdgeKey& edge) -> std::size_t {
        const auto found = midpoint_ids.find(edge);
        if (found == midpoint_ids.end()) {
            throw std::logic_error("A selected adaptive edge has no midpoint.");
        }
        return found->second;
    };
    const auto append = [&](const Triangle2D& child, const std::size_t parent) {
        result.mesh.elements.push_back(child);
        result.parent_elements.push_back(parent);
    };

    for (std::size_t parent = 0; parent < mesh.elements.size(); ++parent) {
        const Triangle2D& t = mesh.elements[parent];
        const std::size_t a = t[0];
        const std::size_t b = t[1];
        const std::size_t c = t[2];
        const EdgeKey ab_key = canonical_edge(a, b);
        const EdgeKey bc_key = canonical_edge(b, c);
        const EdgeKey ca_key = canonical_edge(c, a);
        const bool split_ab = split_edges.contains(ab_key);
        const bool split_bc = split_edges.contains(bc_key);
        const bool split_ca = split_edges.contains(ca_key);
        const unsigned split_count =
            static_cast<unsigned>(split_ab) +
            static_cast<unsigned>(split_bc) +
            static_cast<unsigned>(split_ca);

        if (split_count == 0U) {
            append({a, b, c}, parent);
        } else if (split_count == 1U) {
            if (split_ab) {
                const std::size_t mab = midpoint(ab_key);
                append({a, mab, c}, parent);
                append({mab, b, c}, parent);
            } else if (split_bc) {
                const std::size_t mbc = midpoint(bc_key);
                append({b, mbc, a}, parent);
                append({mbc, c, a}, parent);
            } else {
                const std::size_t mca = midpoint(ca_key);
                append({c, mca, b}, parent);
                append({mca, a, b}, parent);
            }
        } else if (split_count == 2U) {
            if (split_ab && split_bc) {
                const std::size_t mab = midpoint(ab_key);
                const std::size_t mbc = midpoint(bc_key);
                append({a, mab, c}, parent);
                append({mab, b, mbc}, parent);
                append({mab, mbc, c}, parent);
            } else if (split_bc && split_ca) {
                const std::size_t mbc = midpoint(bc_key);
                const std::size_t mca = midpoint(ca_key);
                append({b, mbc, a}, parent);
                append({mbc, c, mca}, parent);
                append({mbc, mca, a}, parent);
            } else {
                const std::size_t mca = midpoint(ca_key);
                const std::size_t mab = midpoint(ab_key);
                append({c, mca, b}, parent);
                append({mca, a, mab}, parent);
                append({mca, mab, b}, parent);
            }
        } else {
            const std::size_t mab = midpoint(ab_key);
            const std::size_t mbc = midpoint(bc_key);
            const std::size_t mca = midpoint(ca_key);
            append({a, mab, mca}, parent);
            append({mab, b, mbc}, parent);
            append({mca, mbc, c}, parent);
            append({mab, mbc, mca}, parent);
        }
    }

    result.boundary_edges.reserve(boundary_edges.size() + split_edges.size());
    for (const BoundaryEdge2D& edge : boundary_edges) {
        const auto found = midpoint_ids.find(canonical_edge(edge[0], edge[1]));
        if (found == midpoint_ids.end()) {
            result.boundary_edges.push_back(edge);
        } else {
            result.boundary_edges.push_back({edge[0], found->second});
            result.boundary_edges.push_back({found->second, edge[1]});
        }
    }

    result.mesh.validate();
    return result;
}


P1ResidualAdaptation2D adapt_p1_diffusion_mesh_by_residual_2d(
    const TriangleMesh2D& mesh,
    const std::vector<BoundaryEdge2D>& boundary_edges,
    const std::span<const double> nodal_values,
    const std::span<const double> element_diffusivity,
    const std::function<double(double, double)>& source,
    const double theta
) {
    DiffusionResidualIndicators2D estimate =
        estimate_p1_residual_indicators_2d(
            mesh, nodal_values, element_diffusivity, source
        );
    std::vector<std::size_t> marked = mark_elements_dorfler_2d(
        estimate.element_indicators, theta
    );
    RefinedTriangleMesh2D refinement =
        refine_triangles_by_marked_elements_2d(
            mesh, boundary_edges, marked
        );
    return {
        std::move(estimate),
        std::move(marked),
        std::move(refinement)
    };
}

} // namespace nova::fem
