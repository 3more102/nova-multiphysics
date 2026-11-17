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
using BoundaryEdge2D = std::array<std::size_t, 2>;

struct TriangleMesh2D {
    std::vector<Point2D> nodes;
    std::vector<Triangle2D> elements;

    // Checks finite coordinates, valid connectivity, and non-degenerate cells.
    void validate() const;
};

struct StructuredRectangleMesh2D {
    TriangleMesh2D mesh;
    std::vector<std::size_t> boundary_nodes;
    std::vector<BoundaryEdge2D> boundary_edges;
};

// Creates a uniform rectangular grid split into two P1 triangles per cell.
// The rectangle is [x_min,x_max] x [y_min,y_max].
[[nodiscard]] StructuredRectangleMesh2D make_rectangle_mesh_2d(
    double x_min,
    double x_max,
    double y_min,
    double y_max,
    std::size_t elements_x,
    std::size_t elements_y
);

} // namespace nova::fem
