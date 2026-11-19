#pragma once

#include "nova/fem/mesh3d.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nova::fem {

// Conforming global 1-to-8 tetrahedral refinement. Every unique coarse edge
// generates one shared midpoint node. Each tetrahedron produces four corner
// children and four children across its central octahedron. The shortest of
// the three octahedral diagonals is chosen deterministically to limit
// unnecessary distortion; no quality improvement is guaranteed.
//
// Uniform refinement is NOT local/adaptive refinement, and it does not
// transfer boundary-condition callbacks, material databases or P2 DOFs.
struct RefinedTetrahedralMesh3D {
    TetrahedralMesh3D mesh;
    std::size_t coarse_node_count{};
    std::size_t coarse_element_count{};
    // New vertex mesh.nodes[coarse_node_count + i] lies on this coarse edge.
    std::vector<std::array<std::size_t, 2>> midpoint_parent_edges;
    // For each child tetrahedron, its parent element index.
    std::vector<std::size_t> parent_elements;
    // Four child faces per supplied exterior face in the supplied face order;
    // each child preserves its parent's outward winding.
    std::vector<BoundaryFace3D> boundary_faces;
    // For each emitted boundary face, its original supplied face index.
    std::vector<std::size_t> parent_boundary_faces;
};

// Supplied boundary_faces may be empty; otherwise they must be distinct
// exterior faces oriented outward from the coarse mesh. The refined mesh
// is revalidated before return. Throws for invalid inputs or overflow.
[[nodiscard]] RefinedTetrahedralMesh3D refine_tetrahedra_uniformly_3d(
    const TetrahedralMesh3D& mesh,
    const std::vector<BoundaryFace3D>& boundary_faces = {}
);

// Exact P1 interpolation at the coarse nodes and new edge midpoints.
// Requires a finite coarse nodal vector with one entry per coarse vertex.
[[nodiscard]] std::vector<double> prolongate_p1_uniform_3d(
    const RefinedTetrahedralMesh3D& refined,
    std::span<const double> coarse_nodal_values
);

// Replicate a scalar element-constant quantity to each of its eight children
// (e.g., per-cell conductivity, heat capacity or source density).
[[nodiscard]] std::vector<double> prolongate_cellwise_constant_uniform_3d(
    const RefinedTetrahedralMesh3D& refined,
    std::span<const double> coarse_cell_values
);

} // namespace nova::fem
