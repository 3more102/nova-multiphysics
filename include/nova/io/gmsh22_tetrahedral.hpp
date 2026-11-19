#pragma once

#include "nova/fem/mesh3d.hpp"

#include <iosfwd>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace nova::io {

// Gmsh MSH 2.2 ASCII: four-node tetrahedra (type 4), exterior triangles
// (type 2), and optional lower-dimensional line/point records.
// Physical IDs are the first element tag, or zero if no tag is present.
// Exterior faces are always returned in NOVA's outward orientation;
// untagged exterior faces are included with physical ID zero.
struct GmshTetrahedralMesh3D {
    fem::TetrahedralMesh3D mesh;
    std::vector<fem::BoundaryFace3D> exterior_faces;
    std::vector<int> exterior_physical_tags;
    std::vector<int> volume_physical_tags;
    // Gmsh $PhysicalNames, keyed by (entity dimension, physical tag).
    // Different dimensions may intentionally use the same tag or name.
    std::map<std::pair<int, int>, std::string> physical_names;
};

// Rejects binary/4.x/high-order or malformed MSH files, missing nodes,
// invalid tetrahedral topology, and triangles that are not exterior faces.
// Node IDs need not be consecutive; output uses consecutive zero-based indices.
[[nodiscard]] GmshTetrahedralMesh3D read_gmsh22_tetrahedra(
    std::istream& input
);

} // namespace nova::io
