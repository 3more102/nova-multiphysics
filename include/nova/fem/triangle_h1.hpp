#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>

namespace nova::fem {

enum class H1TriangleOrder {
    P1,
    P2
};

[[nodiscard]] constexpr std::size_t triangle_h1_dof_count(
    const H1TriangleOrder order
) noexcept {
    return order == H1TriangleOrder::P1 ? 3U : 6U;
}

struct TriangleH1BasisEvaluation {
    std::size_t count{};
    std::array<double, 6> values{};
    std::array<std::array<double, 2>, 6> reference_gradients{};
};

// Evaluates the nodal Lagrange H1 basis on the reference triangle with
// vertices (0,0), (1,0), and (0,1). P1 uses the three vertex nodes.
// P2 adds one node at the midpoint of each reference edge.
//
// Local P2 ordering is:
//   0,1,2 = vertices
//   3     = edge (0,1)
//   4     = edge (1,2)
//   5     = edge (2,0)
[[nodiscard]] inline TriangleH1BasisEvaluation evaluate_triangle_h1_basis(
    const H1TriangleOrder order,
    const double r,
    const double s
) {
    if (order != H1TriangleOrder::P1 && order != H1TriangleOrder::P2) {
        throw std::invalid_argument("Unsupported H1 triangle polynomial order.");
    }

    const double l0 = 1.0 - r - s;
    const double l1 = r;
    const double l2 = s;

    const std::array<std::array<double, 2>, 3> barycentric_gradients{{
        {{-1.0, -1.0}},
        {{ 1.0,  0.0}},
        {{ 0.0,  1.0}}
    }};

    TriangleH1BasisEvaluation result;
    result.count = triangle_h1_dof_count(order);

    if (order == H1TriangleOrder::P1) {
        result.values[0] = l0;
        result.values[1] = l1;
        result.values[2] = l2;
        result.reference_gradients[0] = barycentric_gradients[0];
        result.reference_gradients[1] = barycentric_gradients[1];
        result.reference_gradients[2] = barycentric_gradients[2];
        return result;
    }

    result.values[0] = l0 * (2.0 * l0 - 1.0);
    result.values[1] = l1 * (2.0 * l1 - 1.0);
    result.values[2] = l2 * (2.0 * l2 - 1.0);
    result.values[3] = 4.0 * l0 * l1;
    result.values[4] = 4.0 * l1 * l2;
    result.values[5] = 4.0 * l2 * l0;

    for (std::size_t i = 0; i < 3; ++i) {
        result.reference_gradients[i] = {
            (4.0 * std::array<double, 3>{l0, l1, l2}[i] - 1.0) *
                barycentric_gradients[i][0],
            (4.0 * std::array<double, 3>{l0, l1, l2}[i] - 1.0) *
                barycentric_gradients[i][1]
        };
    }

    result.reference_gradients[3] = {
        4.0 * (l1 * barycentric_gradients[0][0] +
               l0 * barycentric_gradients[1][0]),
        4.0 * (l1 * barycentric_gradients[0][1] +
               l0 * barycentric_gradients[1][1])
    };
    result.reference_gradients[4] = {
        4.0 * (l2 * barycentric_gradients[1][0] +
               l1 * barycentric_gradients[2][0]),
        4.0 * (l2 * barycentric_gradients[1][1] +
               l1 * barycentric_gradients[2][1])
    };
    result.reference_gradients[5] = {
        4.0 * (l0 * barycentric_gradients[2][0] +
               l2 * barycentric_gradients[0][0]),
        4.0 * (l0 * barycentric_gradients[2][1] +
               l2 * barycentric_gradients[0][1])
    };

    return result;
}

} // namespace nova::fem
