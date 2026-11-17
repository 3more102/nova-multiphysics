#include "nova/fem/mesh2d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace nova::fem {

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

    constexpr double degeneracy_factor =
        64.0 * std::numeric_limits<double>::epsilon();

    for (const Triangle2D& element : elements) {
        if (element[0] >= nodes.size() ||
            element[1] >= nodes.size() ||
            element[2] >= nodes.size()) {
            throw std::out_of_range("Triangle connectivity references a missing node.");
        }
        if (element[0] == element[1] ||
            element[1] == element[2] ||
            element[2] == element[0]) {
            throw std::invalid_argument("A triangle cannot repeat a node index.");
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
            (b.x - a.x) * (c.y - a.y) -
            (c.x - a.x) * (b.y - a.y);

        if (!std::isfinite(maximum_edge_squared) || !std::isfinite(determinant)) {
            throw std::invalid_argument("Triangle geometry exceeds supported numeric range.");
        }
        if (maximum_edge == 0.0 ||
            std::abs(determinant) <= degeneracy_factor * maximum_edge_squared) {
            throw std::invalid_argument("Mesh contains a degenerate or near-degenerate triangle.");
        }
    }
}

} // namespace nova::fem
