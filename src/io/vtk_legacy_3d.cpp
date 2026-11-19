#include "nova/io/vtk_legacy_3d.hpp"

#include <cmath>
#include <ios>
#include <limits>
#include <locale>
#include <ostream>
#include <stdexcept>

namespace nova::io {
namespace {

[[nodiscard]] bool valid_field_name(const std::string& name) noexcept {
    if (name.empty()) return false;
    const auto is_letter = [](char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
    };
    if (!is_letter(name.front())) return false;
    for (char c : name) {
        if (!is_letter(c) && (c < '0' || c > '9')) return false;
    }
    return true;
}

// Preserve the caller's locale, precision, and format flags on success or
// failure without buffering an entire (potentially large) tetrahedral mesh.
class StreamFormatGuard {
public:
    explicit StreamFormatGuard(std::ostream& output)
        : output_(output), saved_(nullptr) {
        saved_.copyfmt(output_);
        output_.imbue(std::locale::classic());
        output_.flags(std::ios::dec);
        output_.width(0);
        output_.precision(std::numeric_limits<double>::max_digits10);
    }
    StreamFormatGuard(const StreamFormatGuard&) = delete;
    StreamFormatGuard& operator=(const StreamFormatGuard&) = delete;
    ~StreamFormatGuard() noexcept {
        try { output_.copyfmt(saved_); } catch (...) {}
    }

private:
    std::ostream& output_;
    std::ios saved_;
};

} // namespace

void write_legacy_vtk_scalar_3d(
    std::ostream& output,
    const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& values,
    const std::string& field_name
) {
    if (!output) {
        throw std::runtime_error("VTK output stream is not writable.");
    }
    if (!valid_field_name(field_name)) {
        throw std::invalid_argument("VTK field name must be an ASCII identifier.");
    }
    if (mesh.nodes.size() != values.size()) {
        throw std::invalid_argument("VTK scalar data size must equal mesh vertex count.");
    }
    for (double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("VTK scalar data contains a non-finite value.");
        }
    }
    if (mesh.elements.size() > std::numeric_limits<std::size_t>::max() / 5U) {
        throw std::overflow_error("VTK tetrahedral connectivity size overflow.");
    }
    // Perform all input checks before writing; invalid topology cannot
    // silently yield a plausible but geometrically incorrect visualization.
    mesh.validate();

    StreamFormatGuard guard(output);
    output << "# vtk DataFile Version 3.0\n"
           << "NOVA P1 tetrahedral scalar result\n"
           << "ASCII\n"
           << "DATASET UNSTRUCTURED_GRID\n"
           << "POINTS " << mesh.nodes.size() << " double\n";
    for (const auto& point : mesh.nodes) {
        output << point.x << ' ' << point.y << ' ' << point.z << '\n';
    }

    output << "CELLS " << mesh.elements.size() << ' '
           << (5U * mesh.elements.size()) << '\n';
    for (const auto& cell : mesh.elements) {
        output << "4 " << cell[0] << ' ' << cell[1]
               << ' ' << cell[2] << ' ' << cell[3] << '\n';
    }

    output << "CELL_TYPES " << mesh.elements.size() << '\n';
    for (std::size_t i = 0; i < mesh.elements.size(); ++i) {
        output << "10\n"; // VTK_TETRA (linear P1)
    }

    output << "POINT_DATA " << mesh.nodes.size() << '\n'
           << "SCALARS " << field_name << " double 1\n"
           << "LOOKUP_TABLE default\n";
    for (double value : values) output << value << '\n';

    if (!output) {
        throw std::runtime_error("Failed to write VTK tetrahedral result.");
    }
}

void write_legacy_vtk_conduction_3d(
    std::ostream& output, const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& temperature, std::span<const int> material_ids,
    std::span<const fem::P1TetrahedralFlux3D> cell_fields
) {
    if (material_ids.size() != mesh.elements.size() || cell_fields.size() != mesh.elements.size()) {
        throw std::invalid_argument("VTK conduction requires one material ID and flux per cell.");
    }
    for (const auto& field : cell_fields) {
        if (!std::isfinite(field.flux.x) || !std::isfinite(field.flux.y) || !std::isfinite(field.flux.z)) {
            throw std::invalid_argument("VTK cell heat flux must be finite.");
        }
    }
    write_legacy_vtk_scalar_3d(output, mesh, temperature, "temperature");

    StreamFormatGuard guard(output);
    output << "CELL_DATA " << mesh.elements.size() << '\n'
           << "SCALARS material_id int 1\nLOOKUP_TABLE default\n";
    for (int tag : material_ids) output << tag << '\n';
    output << "VECTORS heat_flux double\n";
    for (const auto& field : cell_fields) {
        output << field.flux.x << ' ' << field.flux.y << ' ' << field.flux.z << '\n';
    }
    if (!output) throw std::runtime_error("Failed to write VTK thermal cell data.");
}

} // namespace nova::io
