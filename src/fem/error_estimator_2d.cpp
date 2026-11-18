#include "nova/fem/error_estimator_2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nova::fem {
namespace {

using EdgeKey = std::pair<std::size_t, std::size_t>;

[[nodiscard]] EdgeKey canonical_edge(
    const std::size_t first,
    const std::size_t second
) {
    return {std::min(first, second), std::max(first, second)};
}

struct ElementGradient {
    double x{};
    double y{};
};

struct TriangleQuadraturePoint {
    double l0{};
    double l1{};
    double l2{};
};

// Symmetric degree-two rule on the reference triangle; weights are 1/3.
constexpr std::array<TriangleQuadraturePoint, 3> degree_two_triangle_rule{{
    {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
    {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}
}};

} // namespace

DiffusionResidualIndicators2D estimate_p1_residual_indicators_2d(
    const TriangleMesh2D& mesh,
    const std::span<const double> nodal_values,
    const std::span<const double> element_diffusivity,
    const std::function<double(double, double)>& source
) {
    mesh.validate();
    if (nodal_values.size() != mesh.nodes.size()) {
        throw std::invalid_argument(
            "P1 residual estimator field size must match the mesh node count."
        );
    }
    if (element_diffusivity.size() != mesh.elements.size()) {
        throw std::invalid_argument(
            "P1 residual estimator diffusivity size must match the element count."
        );
    }
    if (!source) {
        throw std::invalid_argument("P1 residual estimation requires a source function.");
    }

    for (const double value : nodal_values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("P1 residual estimator field must be finite.");
        }
    }
    for (const double value : element_diffusivity) {
        if (!std::isfinite(value) || value <= 0.0) {
            throw std::invalid_argument(
                "P1 residual estimator diffusivity must be finite and positive."
            );
        }
    }

    std::vector<ElementGradient> gradients(mesh.elements.size());
    std::vector<double> indicator_squared(mesh.elements.size(), 0.0);
    std::map<EdgeKey, std::vector<std::size_t>> edge_elements;

    for (std::size_t element_index = 0;
         element_index < mesh.elements.size(); ++element_index) {
        const Triangle2D& element = mesh.elements[element_index];
        const Point2D& p0 = mesh.nodes[element[0]];
        const Point2D& p1 = mesh.nodes[element[1]];
        const Point2D& p2 = mesh.nodes[element[2]];

        const double determinant =
            (p1.x - p0.x) * (p2.y - p0.y) -
            (p2.x - p0.x) * (p1.y - p0.y);
        const double area = 0.5 * std::abs(determinant);
        const std::array<std::array<double, 2>, 3> shape_gradients{{
            {{(p1.y - p2.y) / determinant, (p2.x - p1.x) / determinant}},
            {{(p2.y - p0.y) / determinant, (p0.x - p2.x) / determinant}},
            {{(p0.y - p1.y) / determinant, (p1.x - p0.x) / determinant}}
        }};

        double gradient_x = 0.0;
        double gradient_y = 0.0;
        for (std::size_t local = 0; local < 3U; ++local) {
            gradient_x += nodal_values[element[local]] *
                shape_gradients[local][0];
            gradient_y += nodal_values[element[local]] *
                shape_gradients[local][1];
        }
        if (!std::isfinite(gradient_x) || !std::isfinite(gradient_y) ||
            !std::isfinite(area) || area <= 0.0) {
            throw std::overflow_error("P1 residual estimator element geometry overflowed.");
        }
        gradients[element_index] = {gradient_x, gradient_y};

        const double edge_01 = std::hypot(p1.x - p0.x, p1.y - p0.y);
        const double edge_12 = std::hypot(p2.x - p1.x, p2.y - p1.y);
        const double edge_20 = std::hypot(p0.x - p2.x, p0.y - p2.y);
        const double element_diameter = std::max({edge_01, edge_12, edge_20});
        if (!std::isfinite(element_diameter) || element_diameter <= 0.0) {
            throw std::overflow_error("P1 residual estimator element diameter is invalid.");
        }

        // For P1 and elementwise-constant k, the cell strong residual is f.
        double source_squared_integral = 0.0;
        for (const TriangleQuadraturePoint& point : degree_two_triangle_rule) {
            const double x =
                point.l0 * p0.x + point.l1 * p1.x + point.l2 * p2.x;
            const double y =
                point.l0 * p0.y + point.l1 * p1.y + point.l2 * p2.y;
            const double source_value = source(x, y);
            if (!std::isfinite(source_value)) {
                throw std::invalid_argument(
                    "P1 residual estimator source returned a non-finite value."
                );
            }
            source_squared_integral +=
                (area / 3.0) * source_value * source_value;
        }

        indicator_squared[element_index] =
            element_diameter * element_diameter * source_squared_integral;
        if (!std::isfinite(indicator_squared[element_index])) {
            throw std::overflow_error("P1 cell residual indicator overflowed.");
        }

        edge_elements[canonical_edge(element[0], element[1])].push_back(
            element_index
        );
        edge_elements[canonical_edge(element[1], element[2])].push_back(
            element_index
        );
        edge_elements[canonical_edge(element[2], element[0])].push_back(
            element_index
        );
    }

    for (const auto& [edge, adjacent] : edge_elements) {
        if (adjacent.size() != 2U) {
            continue;
        }
        const Point2D& a = mesh.nodes[edge.first];
        const Point2D& b = mesh.nodes[edge.second];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double edge_length = std::hypot(dx, dy);
        if (!std::isfinite(edge_length) || edge_length <= 0.0) {
            throw std::overflow_error("P1 flux-jump edge length is invalid.");
        }

        // A single canonical normal is used for both cells. The jump sign
        // is immaterial because the edge residual is integrated in squared form.
        const double normal_x = dy / edge_length;
        const double normal_y = -dx / edge_length;
        const std::size_t first = adjacent[0];
        const std::size_t second = adjacent[1];
        const double jump =
            (element_diffusivity[first] * gradients[first].x -
             element_diffusivity[second] * gradients[second].x) * normal_x +
            (element_diffusivity[first] * gradients[first].y -
             element_diffusivity[second] * gradients[second].y) * normal_y;
        const double half_edge_contribution =
            0.5 * edge_length * edge_length * jump * jump;
        if (!std::isfinite(half_edge_contribution)) {
            throw std::overflow_error("P1 flux-jump indicator overflowed.");
        }
        indicator_squared[first] += half_edge_contribution;
        indicator_squared[second] += half_edge_contribution;
        if (!std::isfinite(indicator_squared[first]) ||
            !std::isfinite(indicator_squared[second])) {
            throw std::overflow_error("P1 residual indicator accumulation overflowed.");
        }
    }

    DiffusionResidualIndicators2D result;
    result.element_indicators.reserve(indicator_squared.size());
    double global_squared = 0.0;
    for (const double squared : indicator_squared) {
        result.element_indicators.push_back(std::sqrt(squared));
        global_squared += squared;
        if (!std::isfinite(global_squared)) {
            throw std::overflow_error("Global P1 residual estimate overflowed.");
        }
    }
    result.global_estimate = std::sqrt(global_squared);
    return result;
}

} // namespace nova::fem
