#pragma once

#include "nova/fem/h1_space_2d.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace nova::io {

// Write a scalar P1/P2 finite-element field as legacy ASCII VTK
// UNSTRUCTURED_GRID for viewing in ParaView or another VTK reader.
//
// P1: VTK_TRIANGLE (5), points are the original mesh vertices.
// P2: VTK_QUADRATIC_TRIANGLE (22), points are vertices followed by shared
// edge midpoints in H1DofManager2D order. Both cell types preserve the
// existing element-local nodal ordering.
//
// Throws std::invalid_argument for a mismatched/non-finite field or an invalid
// field name, and std::runtime_error if writing fails.
void write_legacy_vtk_scalar_2d(
    std::ostream& output,
    const fem::H1FunctionSpace2D& space,
    const std::vector<double>& values,
    const std::string& field_name = "temperature"
);

} // namespace nova::io
