#include "nova/io/vtk_legacy_3d.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace nova::io {
namespace {

[[nodiscard]] bool is_ascii_letter(const char c) noexcept {
    return (c >= 'A' && c <= 'Z') ||
           (c >= 'a' && c <= 'z') || c == '_';
}

[[nodiscard]] bool valid_field_name(const std::string& name) noexcept {
    if (name.empty() || !is_ascii_letter(name.front())) {
        return false;
    }
    for (const char c : name) {
        if (!is_ascii_letter(c) && (c < '0' || c > '9')) {
            return false;
        }
    }
    return true;
}

} // namespace

void write_legacy_vtk_scalar_3d(
    std::ostream& output,
    const fem::H1FunctionSpace3D& space,
    const std::vector<double>& values,
    const std::string& field_name
) {
    if (!output) {
        throw std::runtime_error("VTK output stream is not writable.");
    }
    if (!valid_field_name(field_name)) {
        throw std::invalid_argument(
            "VTK scalar field name must be an ASCII identifier."
        );
    }
    if (values.size() != space.dimension()) {
        throw std::invalid_argument(
            "VTK scalar data size must match the 3D H1 space dimension."
        );
    }
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "VTK scalar field contains a non-finite value."
            );
        }
    }

    const auto& mesh = space.mesh();
    const auto& dofs = space.dofs();
    const std::size_t node_count = space.dimension();
    const std::size_t cell_count = dofs.element_count();
    if (mesh.nodes.size() != node_count ||
        mesh.elements.size() != cell_count) {
        throw std::invalid_argument(
            "3D H1 mesh size changed after constructing the function space."
        );
    }
    constexpr std::size_t cell_stride = 5U; // Four vertices + length entry.
    if (cell_count > std::numeric_limits<std::size_t>::max() / cell_stride) {
        throw std::overflow_error("VTK tetrahedral connectivity size overflow.");
    }

    // Construct before writing, so invalid inputs do not leave a partial VTK.
    std::ostringstream vtk;
    vtk.imbue(std::locale::classic());
    vtk << std::setprecision(std::numeric_limits<double>::max_digits10);
    vtk << "# vtk DataFile Version 3.0\n"
        << "NOVA P1 tetrahedral scalar result\n"
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << node_count << " double\n";
    for (const auto& point : mesh.nodes) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z)) {
            throw std::invalid_argument("VTK mesh contains a non-finite coordinate.");
        }
        vtk << point.x << ' ' << point.y << ' ' << point.z << '\n';
    }

    vtk << "CELLS " << cell_count << ' '
        << cell_count * cell_stride << '\n';
    for (std::size_t element = 0; element < cell_count; ++element) {
        const auto& local = dofs.element_dofs(element);
        vtk << "4";
        for (const std::size_t node : local) {
            if (node >= node_count) {
                throw std::invalid_argument("VTK tetrahedron references an invalid node.");
            }
            vtk << ' ' << node;
        }
        vtk << '\n';
    }
    vtk << "CELL_TYPES " << cell_count << '\n';
    for (std::size_t element = 0; element < cell_count; ++element) {
        vtk << "10\n"; // VTK_TETRA
    }

    vtk << "POINT_DATA " << node_count << '\n'
        << "SCALARS " << field_name << " double 1\n"
        << "LOOKUP_TABLE default\n";
    for (const double value : values) {
        vtk << value << '\n';
    }

    if (!vtk) {
        throw std::runtime_error("Failed to construct 3D VTK output.");
    }
    output << vtk.str();
    if (!output) {
        throw std::runtime_error("Failed to write 3D VTK output.");
    }
}

} // namespace nova::io
