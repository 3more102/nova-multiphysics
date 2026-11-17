#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace nova::fem {

struct Point2D {
    double x{};
    double y{};
};

using Triangle2D = std::array<std::size_t, 3>;

struct TriangleMesh2D {
    std::vector<Point2D> nodes;
    std::vector<Triangle2D> elements;

    // Checks finite coordinates, valid connectivity, and non-degenerate cells.
    void validate() const;
};

} // namespace nova::fem
