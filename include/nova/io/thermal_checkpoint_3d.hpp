#pragma once

#include "nova/fem/mesh3d.hpp"

#include <iosfwd>
#include <vector>

namespace nova::io {

// Only the physical restart time and absolute nodal temperatures are stored.
// Material data, boundary conditions, integration scheme and solver settings
// belong to the caller's unchanged TransientHeat3DProblem/Options.
struct TransientHeatCheckpoint3D {
    double time{};
    std::vector<double> temperature;
};

// Versioned, human-readable, lossless-double checkpoint. Stores the mesh's
// exact node coordinates and oriented connectivity so a restart cannot silently
// use a different mesh or DOF ordering. Validates all data before writing.
void write_transient_heat_checkpoint_3d(
    std::ostream& output,
    const fem::TetrahedralMesh3D& mesh,
    double time,
    const std::vector<double>& temperature
);

// Rejects wrong versions, altered meshes, malformed/truncated data,
// non-finite values and trailing tokens. Does not alter the supplied mesh.
[[nodiscard]] TransientHeatCheckpoint3D read_transient_heat_checkpoint_3d(
    std::istream& input,
    const fem::TetrahedralMesh3D& expected_mesh
);

} // namespace nova::io
