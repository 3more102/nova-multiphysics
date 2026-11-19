#pragma once

#include <iosfwd>
#include <string>
#include <vector>

namespace nova::io {

// One legacy VTK file with its physical simulation time.
// Paths are intentionally restricted to safe ASCII basenames, relative
// to the directory containing the generated .pvd collection.
struct VtkTimeSeriesEntry {
    double time{};
    std::string filename;
};

// Write a ParaView XML collection referencing pre-existing legacy .vtk files.
// Requires at least one entry, finite strictly increasing times, unique,
// safe .vtk basenames, and a writable stream. Never owns the referenced files.
void write_pvd_collection(
    std::ostream& output,
    const std::vector<VtkTimeSeriesEntry>& entries
);

} // namespace nova::io
