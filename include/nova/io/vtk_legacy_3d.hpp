#pragma once

#include "nova/fem/mesh3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"

#include <iosfwd>
#include <span>
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

// A continuous nodal temperature field, the original Gmsh volume tags and
// discontinuous physical heat flux vectors in mesh.elements order.
// Optional effective cell heat generation (W/m^3) is written as a scalar
// heat_generation field without changing legacy output when absent.
void write_legacy_vtk_conduction_3d(
    std::ostream& output, const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& temperature, std::span<const int> material_ids,
    std::span<const fem::P1TetrahedralFlux3D> cell_fields,
    std::span<const double> cell_heat_generation = {}
);

// Export a continuous P1 nodal scalar field alongside a discontinuous
// elementwise scalar field, suitable for residual indicators, material
// quantities, or other cell-local diagnostics. The cell field is indexed
// in mesh.elements order. All cell data are validated before any bytes
// are written, and caller stream formatting is preserved.
void write_legacy_vtk_point_and_cell_scalar_3d(
    std::ostream& output,
    const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& nodal_values,
    std::span<const double> cell_values,
    const std::string& point_field_name = "solution",
    const std::string& cell_field_name = "error_indicator"
);

} // namespace nova::io
