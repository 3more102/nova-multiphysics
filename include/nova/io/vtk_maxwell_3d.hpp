#pragma once
#include "nova/physics/maxwell_3d.hpp"
#include <iosfwd>

namespace nova::io {
struct MaxwellCellResult3D {
    physics::MaxwellFieldSample3D centroid;
    double absorbed_power{}; // W
    double mean_absorbed_power_density{}; // W/m^3
};

// E, H and Poynting at centroid equal their cell averages for affine
// lowest-order fields. Mean loss is calculated by exact integration.
[[nodiscard]] std::vector<MaxwellCellResult3D> reconstruct_maxwell_cells_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const physics::MaxwellCellMaterial3D> materials,
    double vacuum_wavenumber, const physics::MaxwellSystem3D& system,
    std::span<const std::complex<double>> edge_field);

// Discontinuous CELL_DATA in SI units. All input validation precedes output.
void write_legacy_vtk_maxwell_3d(
    std::ostream& output, const fem::TetrahedralMesh3D& mesh,
    std::span<const MaxwellCellResult3D> fields);
} // namespace nova::io
