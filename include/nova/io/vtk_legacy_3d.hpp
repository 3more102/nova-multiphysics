#pragma once

#include "nova/fem/h1_space_3d.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace nova::io {

// Export a nodal P1 scalar field on tetrahedra as legacy ASCII VTK
// UNSTRUCTURED_GRID (VTK_TETRA=10) for ParaView and other VTK readers.
// The 3D H1 space has one degree of freedom per mesh vertex.
// Invalid names, mismatched/non-finite values and stream failures throw;
// all input validation occurs before writing to the destination stream.
void write_legacy_vtk_scalar_3d(
    std::ostream& output,
    const fem::H1FunctionSpace3D& space,
    const std::vector<double>& values,
    const std::string& field_name = "temperature"
);

} // namespace nova::io
