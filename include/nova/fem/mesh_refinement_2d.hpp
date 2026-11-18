#pragma once

#include "nova/fem/mesh2d.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nova::fem {

// Mean-ratio shape quality: 4*sqrt(3)*area / (sum of squared edge
// lengths). An equilateral triangle has quality 1; degenerate shapes
// approach 0. Angles are reported in radians.
struct TriangleQuality2D {
    double signed_area{};
    double minimum_angle_radians{};
    double mean_ratio_quality{};
    double longest_edge{};
};

struct TriangleMeshQuality2D {
    std::vector<TriangleQuality2D> elements;
    double minimum_angle_radians{};
    double minimum_mean_ratio_quality{};
    double maximum_edge_length{};
    double total_area{};
};

// Validates connectivity and geometry first; quality is orientation
// independent except for signed_area. No mesh coordinates are changed.
[[nodiscard]] TriangleMeshQuality2D evaluate_mesh_quality_2d(
    const TriangleMesh2D& mesh
);

// Global conforming 1-to-4 refinement: each original unique edge gets
// exactly one midpoint, shared by its adjacent triangles. No hanging
// nodes, independent of element orientation or local vertex ordering.
struct RefinedTriangleMesh2D {
    TriangleMesh2D mesh;
    // Edge (a,b) of the coarse mesh generated new node at
    // mesh.nodes[coarse_node_count + edge_index].
    std::size_t coarse_node_count{};
    std::vector<BoundaryEdge2D> midpoint_parent_edges;
    // Each child triangle maps back to its original parent triangle.
    std::vector<std::size_t> parent_elements;
    // User-supplied exterior edges, split in their original direction.
    std::vector<BoundaryEdge2D> boundary_edges;
};

// boundary_edges can be empty. If provided, every entry must be a
// distinct exterior edge of the coarse mesh.
[[nodiscard]] RefinedTriangleMesh2D refine_triangles_uniformly_2d(
    const TriangleMesh2D& mesh,
    const std::vector<BoundaryEdge2D>& boundary_edges = {}
);

// Interpolate coarse P1 nodal values to the uniformly refined mesh.
// Original vertex values are preserved; midpoint values are averages.
[[nodiscard]] std::vector<double> prolongate_p1_uniform_2d(
    const RefinedTriangleMesh2D& refined,
    std::span<const double> coarse_nodal_values
);


// Interpolate a coarse P2 H1 field onto the uniformly refined P2 space.
// Coarse vertex and edge-midpoint values are preserved. New fine-edge
// midpoint DOFs are evaluated from the parent coarse quadratic element.
[[nodiscard]] std::vector<double> prolongate_p2_uniform_2d(
    const TriangleMesh2D& coarse_mesh,
    const RefinedTriangleMesh2D& refined,
    std::span<const double> coarse_p2_dofs
);

} // namespace nova::fem
