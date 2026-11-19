#pragma once

#include <array>
#include <cstddef>

namespace nova::fem {

// Barycentric shape values; weights include the reference simplex measure.
template<std::size_t Vertices>
struct SimplexQuadraturePoint {
    std::array<double, Vertices> shape;
    double weight;
};

inline constexpr auto tetrahedron_degree_two_quadrature = [] {
    // Permutations of (a,b,b,b), with a=(5+3 sqrt(5))/20, b=(5-sqrt(5))/20.
    constexpr double a = 0.5854101966249685, b = 0.1381966011250105;
    constexpr double w = 1.0 / 24.0;
    return std::array<SimplexQuadraturePoint<4>, 4>{{
        {{a,b,b,b}, w}, {{b,a,b,b}, w}, {{b,b,a,b}, w}, {{b,b,b,a}, w}
    }};
}();

inline constexpr auto triangle_face_degree_two_quadrature = [] {
    constexpr double a = 2.0 / 3.0, b = 1.0 / 6.0, w = 1.0 / 6.0;
    return std::array<SimplexQuadraturePoint<3>, 3>{{
        {{a,b,b}, w}, {{b,a,b}, w}, {{b,b,a}, w}
    }};
}();

} // namespace nova::fem
