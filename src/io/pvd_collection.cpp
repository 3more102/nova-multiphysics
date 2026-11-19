#include "nova/io/pvd_collection.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace nova::io {
namespace {

[[nodiscard]] bool valid_local_filename(const std::string& value) noexcept {
    if (value.empty() || value == "." || value == "..") return false;
    for (unsigned char c : value) {
        const bool letter = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
        const bool digit = c >= '0' && c <= '9';
        if (!letter && !digit && c != '.' && c != '_' && c != '-') return false;
    }
    return true;
}

} // namespace

void write_paraview_collection(
    std::ostream& output, const std::vector<ParaviewTimeFrame>& frames
) {
    if (!output) {
        throw std::runtime_error("ParaView collection stream is not writable.");
    }
    if (frames.empty()) {
        throw std::invalid_argument("ParaView collection requires at least one frame.");
    }
    double previous = -std::numeric_limits<double>::infinity();
    std::unordered_set<std::string> used_filenames;
    for (const auto& frame : frames) {
        if (!std::isfinite(frame.time) || frame.time <= previous ||
            !valid_local_filename(frame.filename) ||
            !used_filenames.insert(frame.filename).second) {
            throw std::invalid_argument("Invalid ParaView frame time or basename.");
        }
        previous = frame.time;
    }

    // Format timestamps separately to preserve the output stream's state and
    // ensure locale-independent XML numbers at round-trip double precision.
    output << "<?xml version=\"1.0\"?>\n"
              "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
              "  <Collection>\n";
    for (const auto& frame : frames) {
        std::ostringstream timestamp;
        timestamp.imbue(std::locale::classic());
        timestamp << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << frame.time;
        output << "    <DataSet timestep=\"" << timestamp.str()
               << "\" group=\"\" part=\"0\" file=\"" << frame.filename
               << "\"/>\n";
    }
    output << "  </Collection>\n</VTKFile>\n";
    if (!output) {
        throw std::runtime_error("Failed to write ParaView time collection.");
    }
}

} // namespace nova::io
