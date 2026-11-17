#include "nova/fem/mesh2d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace nova::fem {
namespace {
using EdgeKey = std::pair<std::size_t, std::size_t>;

[[nodiscard]] EdgeKey canonical_edge(const std::size_t first, const std::size_t second) {
    return {std::min(first, second), std::max(first, second)};
}
} // namespace

void TriangleMesh2D::validate() const {
    if (nodes.empty()) {
        throw std::invalid_argument("A 2D triangle mesh must contain nodes.");
    }
    if (elements.empty()) {
        throw std::invalid_argument("A 2D triangle mesh must contain elements.");
    }
    for (const Point2D& node : nodes) {
        if (!std::isfinite(node.x) || !std::isfinite(node.y)) {
            throw std::invalid_argument("Mesh node coordinates must be finite.");
        }
    }

    constexpr double degeneracy_factor = 64.0 * std::numeric_limits<double>::epsilon();
    std::set<std::array<std::size_t, 3>> unique_elements;

    for (const Triangle2D& element : elements) {
        if (element[0] >= nodes.size() || element[1] >= nodes.size() ||
            element[2] >= nodes.size()) {
            throw std::out_of_range("Triangle connectivity references a missing node.");
        }
        if (element[0] == element[1] || element[1] == element[2] ||
            element[2] == element[0]) {
            throw std::invalid_argument("A triangle cannot repeat a node index.");
        }

        auto canonical_element = element;
        std::sort(canonical_element.begin(), canonical_element.end());
        if (!unique_elements.insert(canonical_element).second) {
            throw std::invalid_argument("Mesh contains a duplicate triangle.");
        }

        const Point2D& a = nodes[element[0]];
        const Point2D& b = nodes[element[1]];
        const Point2D& c = nodes[element[2]];
        const double ab = std::hypot(b.x - a.x, b.y - a.y);
        const double bc = std::hypot(c.x - b.x, c.y - b.y);
        const double ca = std::hypot(a.x - c.x, a.y - c.y);
        const double maximum_edge = std::max({ab, bc, ca});
        const double maximum_edge_squared = maximum_edge * maximum_edge;
        const double determinant =
            (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);

        if (!std::isfinite(maximum_edge_squared) || !std::isfinite(determinant)) {
            throw std::invalid_argument("Triangle geometry exceeds supported numeric range.");
        }
        if (maximum_edge == 0.0 ||
            std::abs(determinant) <= degeneracy_factor * maximum_edge_squared) {
            throw std::invalid_argument("Mesh contains a degenerate or near-degenerate triangle.");
        }
    }

    std::map<EdgeKey, std::size_t> edge_incidence;
    for (const Triangle2D& element : elements) {
        ++edge_incidence[canonical_edge(element[0], element[1])];
        ++edge_incidence[canonical_edge(element[1], element[2])];
        ++edge_incidence[canonical_edge(element[2], element[0])];
    }
    for (const auto& [edge, incidence] : edge_incidence) {
        static_cast<void>(edge);
        if (incidence > 2U) {
            throw std::invalid_argument(
                "Mesh contains a non-manifold edge shared by more than two triangles."
            );
        }
    }
}

StructuredRectangleMesh2D make_rectangle_mesh_2d(
    const double x_min, const double x_max,
    const double y_min, const double y_max,
    const std::size_t elements_x, const std::size_t elements_y
) {
    if (!std::isfinite(x_min) || !std::isfinite(x_max) ||
        !std::isfinite(y_min) || !std::isfinite(y_max) ||
        x_max <= x_min || y_max <= y_min) {
        throw std::invalid_argument(
            "Rectangle bounds must be finite and have positive width and height."
        );
    }
    if (elements_x == 0 || elements_y == 0) {
        throw std::invalid_argument("Rectangle mesh subdivisions must be positive.");
    }

    const std::size_t maximum = std::numeric_limits<std::size_t>::max();
    if (elements_x == maximum || elements_y == maximum) {
        throw std::length_error("Rectangle mesh dimensions exceed supported index range.");
    }
    const std::size_t nodes_x = elements_x + 1;
    const std::size_t nodes_y = elements_y + 1;
    if (nodes_x > maximum / nodes_y ||
        elements_x > maximum / elements_y / 2U ||
        elements_y > maximum / 2U ||
        elements_x > maximum / 2U - elements_y) {
        throw std::length_error("Rectangle mesh size overflows the supported index range.");
    }

    const std::size_t boundary_count = 2U * (elements_x + elements_y);
    StructuredRectangleMesh2D result;
    result.mesh.nodes.reserve(nodes_x * nodes_y);
    result.mesh.elements.reserve(2U * elements_x * elements_y);
    result.boundary_nodes.reserve(boundary_count);
    result.boundary_edges.reserve(boundary_count);

    const double step_x = (x_max - x_min) / static_cast<double>(elements_x);
    const double step_y = (y_max - y_min) / static_cast<double>(elements_y);

    for (std::size_t j = 0; j < nodes_y; ++j) {
        for (std::size_t i = 0; i < nodes_x; ++i) {
            const double x = i == elements_x
                ? x_max : x_min + static_cast<double>(i) * step_x;
            const double y = j == elements_y
                ? y_max : y_min + static_cast<double>(j) * step_y;
            const std::size_t node = j * nodes_x + i;
            result.mesh.nodes.push_back({x, y});
            if (i == 0 || i == elements_x || j == 0 || j == elements_y) {
                result.boundary_nodes.push_back(node);
            }
        }
    }

    for (std::size_t j = 0; j < elements_y; ++j) {
        for (std::size_t i = 0; i < elements_x; ++i) {
            const std::size_t lower_left = j * nodes_x + i;
            const std::size_t lower_right = lower_left + 1;
            const std::size_t upper_left = lower_left + nodes_x;
            const std::size_t upper_right = upper_left + 1;
            result.mesh.elements.push_back({lower_left, lower_right, upper_right});
            result.mesh.elements.push_back({lower_left, upper_right, upper_left});
        }
    }

    // Counter-clockwise boundary edge ordering.
    for (std::size_t i = 0; i < elements_x; ++i) {
        result.boundary_edges.push_back({i, i + 1});
    }
    for (std::size_t j = 0; j < elements_y; ++j) {
        const std::size_t right_bottom = j * nodes_x + elements_x;
        result.boundary_edges.push_back({right_bottom, right_bottom + nodes_x});
    }
    for (std::size_t i = 0; i < elements_x; ++i) {
        const std::size_t top_right = elements_y * nodes_x + elements_x - i;
        result.boundary_edges.push_back({top_right, top_right - 1});
    }
    for (std::size_t j = 0; j < elements_y; ++j) {
        const std::size_t left_top = (elements_y - j) * nodes_x;
        result.boundary_edges.push_back({left_top, left_top - nodes_x});
    }

    result.mesh.validate();
    return result;
}

} // namespace nova::fem
