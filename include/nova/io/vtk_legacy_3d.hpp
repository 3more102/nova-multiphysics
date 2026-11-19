#pragma once

#include "nova/fem/mesh3d.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace nova::io {

// Export a validated P1 tetrahedral mesh and nodal scalar data as legacy
// ASCII VTK UNSTRUCTURED_GRID (VTK_TETRA, cell type 10).
// The scalar values must be indexed by mesh vertex ID, not by cell.
// Rejects invalid field names, non-finite values, invalid meshes and
// unsuccessful output streams before claiming a complete export.
void write_legacy_vtk_scalar_3d(
    std::ostream& output,
    const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& values,
    const std::string& field_name = "temperature"
);

} // namespace nova::io
