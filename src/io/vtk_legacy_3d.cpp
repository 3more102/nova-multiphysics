#include "nova/io/vtk_legacy_3d.hpp"

#include "nova/physics/linear_elasticity_3d.hpp"
#include "nova/physics/package_3d.hpp"

#include <algorithm>
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
        : output_(output), saved_(output.rdbuf()) {
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

void write_scalar(
    std::ostream& output,
    const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& values,
    const std::string& field_name, const physics::LinearElasticity3DResult* mechanics
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
           << "DATASET UNSTRUCTURED_GRID\n";
    if (mechanics) {
        output << "FIELD FieldData 3\n"
               << "elastic_iterations 1 1 double\n" << mechanics->linear_solver.iterations << '\n'
               << "elastic_residual 1 1 double\n" << mechanics->linear_solver.residual_norm << '\n'
               << "elastic_strain_energy_J 1 1 double\n" << mechanics->strain_energy_joules << '\n';
    }
    output << "POINTS " << mesh.nodes.size() << " double\n";
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

void write_tensor(std::ostream& output, const physics::ElasticTensor6& t) {
    // Symmetric Cartesian tensor, tensor shear (xx, yy, zz, xy, yz, xz).
    output << t[0] << ' ' << t[3] << ' ' << t[5] << '\n'
           << t[3] << ' ' << t[1] << ' ' << t[4] << '\n'
           << t[5] << ' ' << t[4] << ' ' << t[2] << '\n';
}

void write_conduction(
    std::ostream& output, const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& temperature, std::span<const int> material_ids,
    std::span<const fem::P1TetrahedralFlux3D> cell_fields,
    std::span<const double> cell_heat_generation, const physics::LinearElasticity3DResult* mechanics,
    const physics::PackageState3D* state=nullptr, const physics::TransportFields3D* transport=nullptr
) {
    if (material_ids.size() != mesh.elements.size() || cell_fields.size() != mesh.elements.size()) {
        throw std::invalid_argument("VTK conduction requires one material ID and flux per cell.");
    }
    if (!cell_heat_generation.empty() && cell_heat_generation.size() != mesh.elements.size()) {
        throw std::invalid_argument("VTK heat generation needs one value per cell.");
    }
    for (double source : cell_heat_generation) {
        if (!std::isfinite(source)) {
            throw std::invalid_argument("VTK heat generation must be finite.");
        }
    }
    for (const auto& field : cell_fields) {
        if (!std::isfinite(field.flux.x) || !std::isfinite(field.flux.y) || !std::isfinite(field.flux.z)) {
            throw std::invalid_argument("VTK cell heat flux must be finite.");
        }
    }
    if (mechanics) {
        const auto finite_field = [](const auto& field) {
            for (const auto& tuple : field) {
                if (!std::all_of(tuple.begin(), tuple.end(), [](double v) { return std::isfinite(v); })) return false;
            }
            return true;
        };
        if (!mechanics->linear_solver.converged || !std::isfinite(mechanics->linear_solver.residual_norm) ||
            !std::isfinite(mechanics->strain_energy_joules) ||
            mechanics->displacement.size() != mesh.nodes.size() ||
            mechanics->cell_stress.size() != mesh.elements.size() ||
            mechanics->cell_strain.size() != mesh.elements.size() ||
            !finite_field(mechanics->displacement) || !finite_field(mechanics->cell_stress) ||
            !finite_field(mechanics->cell_strain)) {
            throw std::invalid_argument("VTK requires a converged, finite mechanical field on the same mesh.");
        }
        for (const auto& field : cell_fields) {
            const auto& g = field.gradient;
            if (!std::isfinite(g.x) || !std::isfinite(g.y) || !std::isfinite(g.z)) {
                throw std::invalid_argument("VTK temperature gradient must be finite.");
            }
        }
    }
    if (transport) {
        const auto valid = [](const auto& values,std::size_t n) {
            return values.size()==n && std::all_of(values.begin(),values.end(),[](double v){return std::isfinite(v);});
        };
        const auto ne=mesh.elements.size(), nn=mesh.nodes.size();
        if (!state || !mechanics || !valid(state->voltage,nn) || transport->electrically_active.size()!=nn ||
            mechanics->support_reaction.size()!=nn || !valid(transport->joule_density,ne) ||
            !valid(transport->conductivity,ne) || !valid(transport->young_modulus,ne) ||
            transport->current_density.size()!=ne) throw std::invalid_argument("Incomplete package VTK fields.");
        for (const auto& j:transport->current_density)
            if (!(std::isfinite(j.x) && std::isfinite(j.y) && std::isfinite(j.z)))
                throw std::invalid_argument("Nonfinite current density.");
        for (const auto& f:mechanics->support_reaction) for (double x:f)
            if (!std::isfinite(x)) throw std::invalid_argument("Nonfinite reaction force.");
    }
    write_scalar(output, mesh, temperature, "temperature", mechanics);

    StreamFormatGuard guard(output);
    if (mechanics) {
        output << "VECTORS displacement double\n";
        for (const auto& u : mechanics->displacement) output << u[0] << ' ' << u[1] << ' ' << u[2] << '\n';
    }
    const auto scalars=[&](const char* name,const auto& values) {
        output << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
        for (auto value:values) output << +value << '\n';
    };
    if (transport) {
        scalars("electric_potential",state->voltage);
        scalars("electrically_active",transport->electrically_active);
        output << "VECTORS support_reaction double\n";
        for (const auto& r:mechanics->support_reaction) output << r[0] << ' ' << r[1] << ' ' << r[2] << '\n';
    }
    output << "CELL_DATA " << mesh.elements.size() << '\n'
           << "SCALARS material_id int 1\nLOOKUP_TABLE default\n";
    for (int tag : material_ids) output << tag << '\n';
    if (!cell_heat_generation.empty()) {
        output << "SCALARS heat_generation double 1\nLOOKUP_TABLE default\n";
        for (double source : cell_heat_generation) output << source << '\n';
    }
    output << "VECTORS heat_flux double\n";
    for (const auto& field : cell_fields) {
        output << field.flux.x << ' ' << field.flux.y << ' ' << field.flux.z << '\n';
    }
    if (mechanics) {
        output << "VECTORS temperature_gradient double\n";
        for (const auto& f : cell_fields) {
            output << f.gradient.x << ' ' << f.gradient.y << ' ' << f.gradient.z << '\n';
        }
        output << "TENSORS strain double\n";
        for (const auto& tensor : mechanics->cell_strain) write_tensor(output, tensor);
        output << "TENSORS mean_stress double\n";
        for (const auto& tensor : mechanics->cell_stress) write_tensor(output, tensor);
    }
    if (transport) {
        scalars("joule_density",transport->joule_density);
        scalars("electrical_conductivity",transport->conductivity);
        scalars("young_modulus",transport->young_modulus);
        output << "VECTORS current_density double\n";
        for (const auto& j:transport->current_density) output << j.x << ' ' << j.y << ' ' << j.z << '\n';
        output << "SCALARS von_mises double 1\nLOOKUP_TABLE default\n";
        for (const auto& s:mechanics->cell_stress) output << physics::von_mises_stress_3d(s) << '\n';
        output << "SCALARS principal_stresses double 3\nLOOKUP_TABLE default\n";
        for (const auto& s:mechanics->cell_stress) {
            const auto v=physics::principal_stresses_3d(s);
            output << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
        }
    }
    if (!output) throw std::runtime_error("Failed to write VTK thermal cell data.");
}

} // namespace

void write_legacy_vtk_scalar_3d(
    std::ostream& output, const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& values, const std::string& field_name
) {
    write_scalar(output, mesh, values, field_name, nullptr);
}

void write_legacy_vtk_conduction_3d(
    std::ostream& output, const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& temperature, std::span<const int> material_ids,
    std::span<const fem::P1TetrahedralFlux3D> cell_fields, std::span<const double> generation
) {
    write_conduction(output, mesh, temperature, material_ids, cell_fields, generation, nullptr);
}

void write_legacy_vtk_thermoelastic_3d(
    std::ostream& output, const fem::TetrahedralMesh3D& mesh,
    const std::vector<double>& temperature, std::span<const int> material_ids,
    std::span<const fem::P1TetrahedralFlux3D> cell_fields,
    const physics::LinearElasticity3DResult& mechanics, std::span<const double> generation
) {
    write_conduction(output, mesh, temperature, material_ids, cell_fields, generation, &mechanics);
}

void write_legacy_vtk_package_3d(std::ostream& out,const fem::TetrahedralMesh3D& mesh,
    std::span<const int> ids,const physics::PackageState3D& state,const physics::PackageStep3D& step) {
    write_conduction(out,mesh,state.temperature,ids,step.transport.heat_flux,{},&step.mechanics,&state,&step.transport);
}

} // namespace nova::io
