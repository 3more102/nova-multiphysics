#pragma once

#include "nova/io/gmsh22_tetrahedral.hpp"

#include <iosfwd>
#include <string_view>

namespace nova::io {

// Import first-order tetrahedral meshes from ASCII Gmsh MSH 2.2 or 4.1.
// MSH 4.1 supports sparse IDs, multiple node/element entity blocks,
// parametric node records, and one physical group per geometric entity.
// Binary files, partitioned/periodic metadata, high-order elements and
// multiple physical tags on one entity are rejected, never silently dropped.
// The resulting mesh and outward exterior faces use NOVA's validated MSH 2.2
// topology and physical-group semantics.
[[nodiscard]] GmshTetrahedralMesh3D read_gmsh_tetrahedra(std::istream& input);

[[nodiscard]] int physical_tag(std::string_view token, int dimension, const GmshTetrahedralMesh3D& mesh);

} // namespace nova::io
