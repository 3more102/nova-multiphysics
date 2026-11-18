#include "nova/io/vtk_legacy_2d.hpp"

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

void write_legacy_vtk_scalar_2d(
    std::ostream& output,
    const fem::H1FunctionSpace2D& space,
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
            "VTK scalar data size must match the H1 space dimension."
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
    const std::size_t nodes_per_cell = space.element_dof_count();
    const std::size_t cell_count = mesh.elements.size();
    const std::size_t cell_stride = nodes_per_cell + 1U;
    if (cell_count > std::numeric_limits<std::size_t>::max() / cell_stride) {
        throw std::overflow_error("VTK connectivity size overflow.");
    }

    // Build the complete record before writing to avoid partially generated
    // output when an input validation error is detected.
    std::ostringstream vtk;
    vtk.imbue(std::locale::classic());
    vtk << std::setprecision(std::numeric_limits<double>::max_digits10);
    vtk << "# vtk DataFile Version 3.0\n"
        << "NOVA H1 scalar result\n"
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << node_count << " double\n";

    for (const auto& point : mesh.nodes) {
        vtk << point.x << ' ' << point.y << " 0\n";
    }
    if (space.order() == fem::H1TriangleOrder::P2) {
        for (const auto& [first, second] : dofs.edge_keys()) {
            const auto& a = mesh.nodes[first];
            const auto& b = mesh.nodes[second];
            vtk << 0.5 * (a.x + b.x) << ' '
                << 0.5 * (a.y + b.y) << " 0\n";
        }
    }

    vtk << "CELLS " << cell_count << ' '
        << cell_count * cell_stride << '\n';
    for (std::size_t element = 0; element < cell_count; ++element) {
        const auto& local = dofs.element_dofs(element);
        vtk << nodes_per_cell;
        for (std::size_t j = 0; j < nodes_per_cell; ++j) {
            vtk << ' ' << local[j];
        }
        vtk << '\n';
    }

    const int vtk_cell_type =
        space.order() == fem::H1TriangleOrder::P1 ? 5 : 22;
    vtk << "CELL_TYPES " << cell_count << '\n';
    for (std::size_t element = 0; element < cell_count; ++element) {
        vtk << vtk_cell_type << '\n';
    }

    vtk << "POINT_DATA " << node_count << '\n'
        << "SCALARS " << field_name << " double 1\n"
        << "LOOKUP_TABLE default\n";
    for (const double value : values) {
        vtk << value << '\n';
    }

    if (!vtk) {
        throw std::runtime_error("Failed to construct VTK output.");
    }
    output << vtk.str();
    if (!output) {
        throw std::runtime_error("Failed to write VTK output.");
    }
}

} // namespace nova::io
