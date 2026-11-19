#include "nova/io/vtk_pvd.hpp"

#include <cmath>
#include <iomanip>
#include <ios>
#include <limits>
#include <locale>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace nova::io {
namespace {

void validate_filename(std::string_view filename) {
    if (filename.empty() || filename.front() == '/' ||
        filename.back() == '/' ||
        !(filename.ends_with(".vtk") || filename.ends_with(".vtu"))) {
        throw std::invalid_argument("PVD requires a relative .vtk or .vtu file path.");
    }

    for (const unsigned char c : filename) {
        if (c < 0x20U || c >= 0x7fU ||
            c == '\\' || c == ':' || c == '<' || c == '>' ||
            c == '"' || c == '|' || c == '?' || c == '*') {
            throw std::invalid_argument("PVD snapshot path contains a nonportable character.");
        }
    }

    std::size_t start = 0;
    while (start < filename.size()) {
        const auto end = filename.find('/', start);
        const auto segment = filename.substr(
            start, end == std::string_view::npos ? end : end - start);
        if (segment.empty() || segment == "." || segment == "..") {
            throw std::invalid_argument("PVD snapshot path has an empty or unsafe segment.");
        }
        if (end == std::string_view::npos) break;
        start = end + 1;
    }
}

std::string escape_attribute(std::string_view value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char c : value) {
        switch (c) {
        case '&': escaped += "&amp;"; break;
        case '<': escaped += "&lt;"; break;
        case '>': escaped += "&gt;"; break;
        case '"': escaped += "&quot;"; break;
        case '\'': escaped += "&apos;"; break;
        default: escaped.push_back(c); break;
        }
    }
    return escaped;
}

} // namespace

void write_vtk_pvd_collection(
    std::ostream& output,
    std::span<const VtkTimeSeriesEntry> entries
) {
    if (!output) {
        throw std::runtime_error("PVD output stream is not writable.");
    }
    if (entries.empty()) {
        throw std::invalid_argument("PVD collection needs at least one timestep.");
    }
    // Validate the complete series before emitting even the XML declaration.
    std::set<std::string> filenames;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto& entry = entries[i];
        if (!std::isfinite(entry.time)) {
            throw std::invalid_argument("PVD timestep must be finite.");
        }
        if (i > 0 && !(entry.time > entries[i - 1].time)) {
            throw std::invalid_argument("PVD timesteps must strictly increase.");
        }
        validate_filename(entry.filename);
        if (!filenames.insert(entry.filename).second) {
            throw std::invalid_argument("PVD snapshot filenames must be unique.");
        }
    }

    // Buffer only the small collection index, not its referenced VTK data.
    // The classic locale guarantees a dot decimal separator on every platform.
    std::ostringstream xml;
    xml.imbue(std::locale::classic());
    xml << std::setprecision(std::numeric_limits<double>::max_digits10);
    xml << "<?xml version=\"1.0\"?>\n"
           "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           "  <Collection>\n";
    for (const auto& entry : entries) {
        xml << "    <DataSet timestep=\"" << entry.time
            << "\" group=\"\" part=\"0\" file=\""
            << escape_attribute(entry.filename) << "\"/>\n";
    }
    xml << "  </Collection>\n</VTKFile>\n";
    if (!xml) throw std::runtime_error("Failed to format the PVD collection.");
    output << xml.str();
    if (!output) throw std::runtime_error("Failed to write PVD collection.");
}

} // namespace nova::io
