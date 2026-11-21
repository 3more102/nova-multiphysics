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

// Deterministic Dörfler bulk marking. Indicators must be finite and
// nonnegative; theta must lie in (0,1]. Returns the smallest descending-
// indicator prefix whose squared sum reaches theta times the global sum.
// Ties are resolved by ascending element index; an all-zero field marks none.
[[nodiscard]] std::vector<std::size_t> mark_dorfler_elements_3d(
    std::span<const double> element_indicators,
    double theta
);

// One-level conforming local refinement. Each marked tetrahedron contributes
// its six edges to a global split set; every incident tetrahedron is bisected
// along those shared edges in deterministic order. This subdivides shared
// faces identically and introduces no hanging faces. It is not a multilevel
// grading or mesh-quality optimization algorithm.
[[nodiscard]] RefinedTetrahedralMesh3D refine_tetrahedra_adaptively_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const std::size_t> marked_elements,
    const std::vector<BoundaryFace3D>& boundary_faces = {}
);

// One indicator-driven adaptation cycle: validate indicator count, apply
// deterministic Dörfler marking, and perform conforming local refinement.
// An all-zero indicator field produces an unchanged mesh. The returned
// marked indices preserve the marking priority order for diagnostics.
// This does not estimate errors, transfer physics models or resolve a PDE.
struct IndicatorRefinementResult3D {
    std::vector<std::size_t> marked_elements;
    RefinedTetrahedralMesh3D refinement;
};

[[nodiscard]] IndicatorRefinementResult3D refine_tetrahedra_by_indicators_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const double> element_indicators,
    double theta,
    const std::vector<BoundaryFace3D>& boundary_faces = {}
);

// Exact P1 interpolation at the coarse nodes and new edge midpoints for any
// RefinedTetrahedralMesh3D result (uniform or local).
[[nodiscard]] std::vector<double> prolongate_p1_refined_3d(
    const RefinedTetrahedralMesh3D& refined,
    std::span<const double> coarse_nodal_values
);

// Replicate a scalar element-constant quantity to every child using the
// parent-element map (uniform or local refinement).
[[nodiscard]] std::vector<double> prolongate_cellwise_constant_refined_3d(
    const RefinedTetrahedralMesh3D& refined,
    std::span<const double> coarse_cell_values
);

// Backward-compatible uniform-refinement transfer API.
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
