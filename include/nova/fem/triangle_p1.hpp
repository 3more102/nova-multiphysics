#pragma once

#include <array>

namespace nova::fem {

// Six-point positive degree-four quadrature on the reference triangle
// with vertices (0,0), (1,0), and (0,1). Weights sum to 1/2.
struct TriangleQuadraturePoint {
    double r{};
    double s{};
    double weight{};
};

inline constexpr double triangle_q4_a = 0.445948490915965;
inline constexpr double triangle_q4_b = 0.108103018168070;
inline constexpr double triangle_q4_c = 0.091576213509771;
inline constexpr double triangle_q4_d = 0.816847572980459;
inline constexpr double triangle_q4_w1 = 0.111690794839005;
inline constexpr double triangle_q4_w2 = 0.054975871827661;

inline constexpr std::array<TriangleQuadraturePoint, 6>
triangle_degree_four_quadrature{{
    {triangle_q4_a, triangle_q4_a, triangle_q4_w1},
    {triangle_q4_b, triangle_q4_a, triangle_q4_w1},
    {triangle_q4_a, triangle_q4_b, triangle_q4_w1},
    {triangle_q4_c, triangle_q4_c, triangle_q4_w2},
    {triangle_q4_d, triangle_q4_c, triangle_q4_w2},
    {triangle_q4_c, triangle_q4_d, triangle_q4_w2}
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
