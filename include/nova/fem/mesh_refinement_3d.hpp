#pragma once

#include "nova/fem/mesh3d.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace nova::fem {

// Uniform conforming 1-to-8 refinement of straight-sided P1 tetrahedra.
// Original vertices keep their indices; every unique original edge receives
// exactly one new midpoint, even when shared by many tetrahedra.
struct UniformTetrahedralRefinement3D {
    TetrahedralMesh3D mesh;
    std::size_t original_node_count{};
    // Entry i belongs to node original_node_count + i.
    std::vector<std::array<std::size_t, 2>> midpoint_parent_edges;
    // Entry i identifies the original tetrahedron of refined element i.
    std::vector<std::size_t> parent_elements;
    // Four outward child faces for every original exterior face.
    std::vector<BoundaryFace3D> boundary_faces;
    // Index into extract_boundary_faces_3d(original_mesh) for each child face.
    // This mapping permits retaining exterior surface tags and boundary data.
    std::vector<std::size_t> parent_boundary_faces;
};

// Splits each input tetrahedron into 4 corner tetrahedra and 4 tetrahedra
// across an interior octahedron diagonal. The face subdivision is independent
// of neighboring cell ordering, so refinement is conforming across shared
// coarse faces. The input and refined meshes are both validated.
[[nodiscard]] UniformTetrahedralRefinement3D refine_tetrahedra_uniform_3d(
    const TetrahedralMesh3D& mesh
);

// Exact nodal P1 interpolation: original DOFs are copied; shared new edge
// midpoint DOFs are the mean of their coarse endpoint values.
// Throws on wrong field sizes, nonfinite data or damaged refinement metadata.
[[nodiscard]] std::vector<double> prolongate_tetrahedral_p1_3d(
    const UniformTetrahedralRefinement3D& refinement,
    const std::vector<double>& coarse_values
);

} // namespace nova::fem
