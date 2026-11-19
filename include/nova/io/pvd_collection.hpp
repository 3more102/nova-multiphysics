#pragma once

#include <iosfwd>
#include <string>
#include <vector>

namespace nova::io {

struct ParaviewTimeFrame {
    double time{};
    // Basename only, containing ASCII letters, digits, '.', '_' or '-'.
    // The caller must write this VTK file alongside the PVD collection.
    std::string filename;
};

// Write a ParaView XML collection that references independently written
// VTK datasets. Requires nonempty, strictly increasing finite times and
// safe, local ASCII frame filenames. Validates before writing any bytes.
void write_paraview_collection(
    std::ostream& output,
    const std::vector<ParaviewTimeFrame>& frames
);

} // namespace nova::io
