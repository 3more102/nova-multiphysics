#pragma once

#include <array>

namespace nova::fem {

// Three-point degree-two quadrature on the reference triangle
// with vertices (0,0), (1,0), and (0,1). Weights sum to 1/2.
struct TriangleQuadraturePoint {
    double r{};
    double s{};
    double weight{};
};

inline constexpr std::array<TriangleQuadraturePoint, 3>
triangle_degree_two_quadrature{{
    {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
    {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}
}};

// Linear P1 shape functions in barycentric coordinates.
[[nodiscard]] constexpr std::array<double, 3>
triangle_p1_shape_values(const double r, const double s) noexcept {
    return {1.0 - r - s, r, s};
}

inline constexpr std::array<std::array<double, 2>, 3>
triangle_p1_reference_gradients{{
    {{-1.0, -1.0}},
    {{ 1.0,  0.0}},
    {{ 0.0,  1.0}}
}};

} // namespace nova::fem
