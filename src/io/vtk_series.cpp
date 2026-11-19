#include "nova/io/vtk_series.hpp"

#include <cmath>
#include <iomanip>
#include <ios>
#include <limits>
#include <locale>
#include <ostream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace nova::io {
namespace {

[[nodiscard]] bool is_safe_vtk_basename(const std::string& filename) {
    if (filename.size() < 5 || filename.front() == '.' ||
        filename.find("..") != std::string::npos ||
        filename.compare(filename.size() - 4, 4, ".vtk") != 0) {
        return false;
    }
    for (const char c : filename) {
        const bool alpha = (c >= 'a' && c <= 'z') ||
                           (c >= 'A' && c <= 'Z');
        const bool digit = c >= '0' && c <= '9';
        if (!alpha && !digit && c != '_' && c != '-' && c != '.') {
            return false;
        }
    }
    return true;
}

class StreamFormatGuard {
public:
    explicit StreamFormatGuard(std::ostream& stream)
        : stream_(stream), saved_(nullptr) {
        saved_.copyfmt(stream);
    }
    StreamFormatGuard(const StreamFormatGuard&) = delete;
    StreamFormatGuard& operator=(const StreamFormatGuard&) = delete;
    ~StreamFormatGuard() noexcept {
        try { stream_.copyfmt(saved_); } catch (...) {}
    }

private:
    std::ostream& stream_;
    std::ios saved_;
};

} // namespace

void write_pvd_collection(
    std::ostream& output,
    const std::vector<VtkTimeSeriesEntry>& entries
) {
    if (!output) {
        throw std::runtime_error("PVD collection stream is not writable.");
    }
    if (entries.empty()) {
        throw std::invalid_argument("PVD collection requires at least one timestep.");
    }
    std::unordered_set<std::string> seen;
    double previous_time = 0.0;
    bool first = true;
    // Validate every entry before writing a partial collection.
    for (const auto& entry : entries) {
        if (!std::isfinite(entry.time) ||
            (!first && !(entry.time > previous_time))) {
            throw std::invalid_argument("PVD timesteps must be finite and strictly increasing.");
        }
        if (!is_safe_vtk_basename(entry.filename) ||
            !seen.insert(entry.filename).second) {
            throw std::invalid_argument("PVD requires unique, safe .vtk basenames.");
        }
        previous_time = entry.time;
        first = false;
    }

    StreamFormatGuard guard(output);
    output.imbue(std::locale::classic());
    output.unsetf(std::ios::floatfield); // Do not inherit fixed/scientific from the caller.
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           << "  <Collection>\n";
    for (const auto& entry : entries) {
        output << "    <DataSet timestep=\"" << entry.time
               << "\" group=\"\" part=\"0\" file=\"" << entry.filename
               << "\"/>\n";
    }
    output << "  </Collection>\n</VTKFile>\n";
    if (!output) {
        throw std::runtime_error("Failed to write PVD time-series collection.");
    }
}

} // namespace nova::io
