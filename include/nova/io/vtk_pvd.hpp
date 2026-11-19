#pragma once

#include <iosfwd>
#include <span>
#include <string>

namespace nova::io {

// Each entry refers to a separately written legacy .vtk or XML .vtu dataset.
// Paths are relative to the directory containing the .pvd collection.
struct VtkTimeSeriesEntry {
    double time{};
    std::string filename;
};

// Write a ParaView/VTK XML Collection (.pvd) manifest for an existing sequence
// of VTK snapshots. Entries must have distinct relative filenames and strictly
// increasing finite physical times. Filenames are validated for portable ASCII
// paths, traversal is rejected, and XML attribute contents are escaped.
// Prevalidates every entry before writing; throws on bad or failed streams.
// This function writes the manifest only, not the referenced VTK snapshots.
void write_vtk_pvd_collection(
    std::ostream& output,
    std::span<const VtkTimeSeriesEntry> entries
);

} // namespace nova::io
