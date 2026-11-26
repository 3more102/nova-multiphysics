#include "nova/io/vtk_maxwell_3d.hpp"
#include <algorithm>
#include <cmath>
#include <ios>
#include <limits>
#include <locale>
#include <ostream>
#include <stdexcept>

namespace nova::io {
namespace {
class FormatGuard {
public:
    explicit FormatGuard(std::ostream& out): out_(out), saved_(nullptr) {
        saved_.copyfmt(out_);
        out_.imbue(std::locale::classic());
        out_.flags(std::ios::dec);
        out_.width(0);
        out_.precision(std::numeric_limits<double>::max_digits10);
    }
    ~FormatGuard() noexcept { try { out_.copyfmt(saved_); } catch (...) {} }
    FormatGuard(const FormatGuard&)=delete;
    FormatGuard& operator=(const FormatGuard&)=delete;
private:
    std::ostream& out_;
    std::ios saved_;
};
void validate_field(const MaxwellCellResult3D& f) {
    for (const auto& v : {f.centroid.electric_field,f.centroid.magnetic_field})
        for (const auto z : v)
            if (!std::isfinite(z.real()) || !std::isfinite(z.imag()))
                throw std::invalid_argument("Nonfinite Maxwell field.");
    for (double v : f.centroid.time_averaged_poynting)
        if (!std::isfinite(v)) throw std::invalid_argument("Nonfinite Poynting flux.");
    for (double v : {f.centroid.absorbed_power_density,f.absorbed_power,
                     f.mean_absorbed_power_density})
        if (!(v>=0.0) || !std::isfinite(v))
            throw std::invalid_argument("Invalid passive Maxwell loss.");
}
} // namespace

std::vector<MaxwellCellResult3D> reconstruct_maxwell_cells_3d(
    const fem::TetrahedralMesh3D& mesh,
    std::span<const physics::MaxwellCellMaterial3D> materials,
    double vacuum_wavenumber, const physics::MaxwellSystem3D& system,
    std::span<const std::complex<double>> edge_field) {
    // Existing exact integration checks mesh and entire edge ordering.
    const auto powers=physics::integrate_maxwell_absorption_3d(
        mesh,materials,vacuum_wavenumber,system,edge_field);
    std::vector<MaxwellCellResult3D> fields;
    fields.reserve(mesh.elements.size());
    for (std::size_t c=0;c<mesh.elements.size();++c) {
        const auto cell=fem::make_nedelec_tetrahedron_3d(mesh,c);
        std::array<std::complex<double>,6> local{};
        for (std::size_t e=0;e<local.size();++e) {
            const auto it=std::lower_bound(system.edges.begin(),system.edges.end(),cell.edges[e]);
            if (it==system.edges.end() || *it!=cell.edges[e])
                throw std::invalid_argument("Missing Maxwell edge.");
            local[e]=edge_field[static_cast<std::size_t>(it-system.edges.begin())];
        }
        MaxwellCellResult3D f{
            physics::evaluate_maxwell_field_3d(cell,local,materials[c],vacuum_wavenumber,
                                              0.25,0.25,0.25),
            powers[c],powers[c]/cell.volume};
        validate_field(f);
        fields.push_back(f);
    }
    return fields;
}

void write_legacy_vtk_maxwell_3d(std::ostream& output,
    const fem::TetrahedralMesh3D& mesh, std::span<const MaxwellCellResult3D> fields) {
    if (!output) throw std::runtime_error("Maxwell VTK stream not writable.");
    if (fields.size()!=mesh.elements.size())
        throw std::invalid_argument("Maxwell VTK field count mismatch.");
    constexpr auto limit=static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (mesh.nodes.size()>limit || mesh.elements.size()>limit/5)
        throw std::overflow_error("Maxwell VTK connectivity exceeds legacy limits.");
    mesh.validate();
    for (const auto& f : fields) validate_field(f);
    FormatGuard guard(output);
    output << "# vtk DataFile Version 3.0\nNOVA Maxwell cell fields (SI)\nASCII\n"
           << "DATASET UNSTRUCTURED_GRID\nPOINTS " << mesh.nodes.size() << " double\n";
    for (const auto& p : mesh.nodes) output<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    output<<"CELLS "<<mesh.elements.size()<<' '<<5*mesh.elements.size()<<'\n';
    for (const auto& cell : mesh.elements)
        output<<"4 "<<cell[0]<<' '<<cell[1]<<' '<<cell[2]<<' '<<cell[3]<<'\n';
    output<<"CELL_TYPES "<<mesh.elements.size()<<'\n';
    for (std::size_t c=0;c<mesh.elements.size();++c) output<<"10\n";
    output<<"CELL_DATA "<<fields.size()<<'\n';
    const auto complex_vector=[&](const char* name, bool electric, bool imaginary) {
        output<<"VECTORS "<<name<<" double\n";
        for (const auto& f : fields) {
            const auto& v=electric?f.centroid.electric_field:f.centroid.magnetic_field;
            for (std::size_t d=0;d<3;++d)
                output<<(imaginary?v[d].imag():v[d].real())<<(d==2?'\n':' ');
        }
    };
    complex_vector("electric_real",true,false);
    complex_vector("electric_imag",true,true);
    complex_vector("magnetic_real",false,false);
    complex_vector("magnetic_imag",false,true);
    output<<"VECTORS poynting_mean double\n";
    for (const auto& f : fields) {
        const auto& v=f.centroid.time_averaged_poynting;
        output<<v[0]<<' '<<v[1]<<' '<<v[2]<<'\n';
    }
    output<<"SCALARS absorbed_power double 1\nLOOKUP_TABLE default\n";
    for (const auto& f : fields) output<<f.absorbed_power<<'\n';
    output<<"SCALARS mean_absorbed_power_density double 1\nLOOKUP_TABLE default\n";
    for (const auto& f : fields) output<<f.mean_absorbed_power_density<<'\n';
    if (!output) throw std::runtime_error("Failed to write Maxwell VTK.");
}
} // namespace nova::io
