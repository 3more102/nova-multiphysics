#include "nova/io/thermal_project.hpp"

#include "atomic_file.hpp"
#include "nova/io/vtk_legacy_3d.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace nova::io {
namespace {

using namespace nova::model;
constexpr std::array<std::string_view, 6> side_names{"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
using Preconditioner = linalg::PreconditionerType;

std::string_view preconditioner_name(Preconditioner value) {
    switch (value) {
    case Preconditioner::None: return "none";
    case Preconditioner::Jacobi: return "jacobi";
    case Preconditioner::IncompleteCholesky0: return "ic0";
    }
    throw std::invalid_argument("Unknown CG preconditioner.");
}

std::string bounded_text(std::istream& input) {
    std::string text(max_thermal_project_bytes + 1, '\0');
    input.read(text.data(), static_cast<std::streamsize>(text.size()));
    if (input.bad() || (input.fail() && !input.eof())) {
        throw std::runtime_error("Failed to read thermal project.");
    }
    text.resize(static_cast<std::size_t>(input.gcount()));
    if (text.size() > max_thermal_project_bytes) throw std::invalid_argument("Thermal project exceeds 64 KiB.");
    for (unsigned char ch : text) {
        if (ch != '\n' && ch != '\r' && ch != '\t' && (ch < 32 || ch > 126)) {
            throw std::invalid_argument("Thermal project contains a non-ASCII or control byte.");
        }
    }
    return text;
}

struct Record {
    std::size_t line;
    std::vector<std::string> words;

    [[noreturn]] void fail(const std::string& message) const {
        throw std::invalid_argument("Thermal project line " + std::to_string(line) + ": " + message);
    }
    void count(std::size_t n) const {
        if (words.size() != n) fail("Incorrect field count for '" + words.front() + "'.");
    }
    void expect(std::size_t i, std::string_view value) const {
        if (words.at(i) != value) fail("Expected '" + std::string(value) + "', got '" + words[i] + "'.");
    }
    template<class T> T number(std::size_t i) const {
        const auto& word = words.at(i);
        T value{};
        const auto end = word.data() + word.size();
        const auto result = std::from_chars(word.data(), end, value);
        if (result.ec != std::errc{} || result.ptr != end) fail("Invalid number '" + word + "'.");
        return value;
    }
};

class ProjectReader {
public:
    explicit ProjectReader(std::istream& input) : source_(bounded_text(input)) {}

    std::optional<Record> next() {
        for (std::string text; std::getline(source_, text);) {
            ++line_;
            const auto comment = text.find('#');
            if (comment != std::string::npos) text.resize(comment);
            std::istringstream fields(text);
            fields.imbue(std::locale::classic());
            Record row{line_, {}};
            for (std::string word; fields >> word;) row.words.push_back(std::move(word));
            if (!row.words.empty()) return row;
        }
        return std::nullopt;
    }
    Record take() {
        if (auto next_row = next()) return std::move(*next_row);
        throw std::invalid_argument("Unexpected end of thermal project after line " + std::to_string(line_) + ".");
    }
    Record take(std::string_view keyword, std::size_t length) {
        auto row = take();
        row.expect(0, keyword);
        row.count(length);
        return row;
    }

private:
    std::istringstream source_;
    std::size_t line_{};
};

void read_boundary(const Record& row, ThermalProject& project) {
    row.expect(0, "boundary");
    if (row.words.size() < 6) row.fail("Incomplete boundary.");
    row.expect(1, project.study_id);
    row.expect(2, project.geometry.id);
    const auto side = std::find(side_names.begin(), side_names.end(), row.words[3]);
    if (side == side_names.end()) row.fail("Unknown box side '" + row.words[3] + "'.");
    ThermalBoundary condition;
    if (row.words[4] == "temperature") {
        row.count(6);
        condition = FixedTemperature{row.number<double>(5)};
    } else if (row.words[4] == "outward_flux") {
        row.count(6);
        condition = OutwardHeatFlux{row.number<double>(5)};
    } else if (row.words[4] == "convection") {
        row.count(7);
        condition = Convection{row.number<double>(5), row.number<double>(6)};
    } else {
        row.fail("Unknown thermal boundary type '" + row.words[4] + "'.");
    }
    if (!project.boundaries.emplace(static_cast<BoxSide>(side - side_names.begin()), condition).second) {
        row.fail("Duplicate assignment on side '" + row.words[3] + "'.");
    }
}

std::string canonical(const ThermalProject& p) {
    validate_thermal_project(p);
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "NOVA_PROJECT " << (p.material.conductivity_tensor_w_m_k ? 2 : 1)
        << "\\nproject " << p.id << "\\ngeometry " << p.geometry.id << " box";
    for (double length : p.geometry.lengths_m) out << ' ' << length;
    out << "\nmesh " << p.mesh.id << ' ' << p.geometry.id << " structured_tet";
    for (auto count : p.mesh.divisions) out << ' ' << count;
    out << "\nmaterial " << p.material.id << ' ';
    if (p.material.conductivity_tensor_w_m_k) {
        const auto& k = *p.material.conductivity_tensor_w_m_k;
        out << "anisotropic " << k.xx << ' ' << k.yy << ' ' << k.zz << ' '
            << k.xy << ' ' << k.xz << ' ' << k.yz;
    } else {
        out << "isotropic " << p.material.conductivity_w_m_k;
    }
    out << "\nstudy " << p.study_id << " steady_heat " << p.mesh.id << ' ' << p.material.id
        << "\nsource " << p.study_id << ' ' << p.generation_w_m3
        << "\nsolver " << p.study_id << " cg " << preconditioner_name(p.solver.preconditioner)
        << ' ' << p.solver.max_iterations << ' ' << p.solver.relative_tolerance << ' '
        << p.solver.absolute_tolerance << '\n';
    for (const auto& [side, condition] : p.boundaries) {
        out << "boundary " << p.study_id << ' ' << p.geometry.id << ' '
            << side_names.at(static_cast<std::size_t>(side)) << ' ';
        if (const auto* fixed = std::get_if<FixedTemperature>(&condition)) {
            out << "temperature " << fixed->kelvin;
        } else if (const auto* flux = std::get_if<OutwardHeatFlux>(&condition)) {
            out << "outward_flux " << flux->watts_per_m2;
        } else {
            const auto& film = std::get<Convection>(condition);
            out << "convection " << film.film_w_m2_k << ' ' << film.ambient_kelvin;
        }
        out << '\n';
    }
    out << "end\n";
    if (!out) throw std::runtime_error("Failed to serialize thermal project.");
    return out.str();
}

} // namespace

model::ThermalProject read_thermal_project(std::istream& input) {
    ProjectReader reader(input);
    const auto header = reader.take("NOVA_PROJECT", 2);
    const auto& schema = header.words[1];
    if (schema != "1" && schema != "2") {
        header.fail("Unsupported NOVA project schema '" + schema + "'.");
    }

    model::ThermalProject p;
    p.id = reader.take("project", 2).words[1];
    const auto geometry = reader.take("geometry", 6);
    geometry.expect(2, "box");
    p.geometry.id = geometry.words[1];
    const auto mesh = reader.take("mesh", 7);
    mesh.expect(2, p.geometry.id);
    mesh.expect(3, "structured_tet");
    p.mesh.id = mesh.words[1];
    for (std::size_t axis = 0; axis < 3; ++axis) {
        p.geometry.lengths_m[axis] = geometry.number<double>(3 + axis);
        p.mesh.divisions[axis] = mesh.number<std::size_t>(4 + axis);
    }
    const auto material = reader.take();
    material.expect(0, "material");
    if (material.words.size() < 3) material.fail("Incomplete material record.");
    p.material.id = material.words[1];
    if (material.words[2] == "isotropic") {
        material.count(4);
        p.material.conductivity_w_m_k = material.number<double>(3);
        p.material.conductivity_tensor_w_m_k.reset();
    } else if (material.words[2] == "anisotropic") {
        if (schema != "2") material.fail("Anisotropic conductivity requires NOVA_PROJECT 2.");
        material.count(9);
        p.material.conductivity_tensor_w_m_k = fem::SymmetricTensor3D{
            material.number<double>(3), material.number<double>(4), material.number<double>(5),
            material.number<double>(6), material.number<double>(7), material.number<double>(8)};
    } else {
        material.fail("Unknown thermal material type '" + material.words[2] + "'.");
    }
    const auto study = reader.take("study", 5);
    study.expect(2, "steady_heat");
    study.expect(3, p.mesh.id);
    study.expect(4, p.material.id);
    p.study_id = study.words[1];
    const auto source = reader.take("source", 3);
    source.expect(1, p.study_id);
    p.generation_w_m3 = source.number<double>(2);
    const auto solver = reader.take("solver", 7);
    solver.expect(1, p.study_id);
    solver.expect(2, "cg");
    const auto& name = solver.words[3];
    if (name == "jacobi") p.solver.preconditioner = Preconditioner::Jacobi;
    else if (name == "ic0") p.solver.preconditioner = Preconditioner::IncompleteCholesky0;
    else if (name == "none") p.solver.preconditioner = Preconditioner::None;
    else solver.fail("Unknown CG preconditioner '" + name + "'.");
    p.solver.max_iterations = solver.number<std::size_t>(4);
    p.solver.relative_tolerance = solver.number<double>(5);
    p.solver.absolute_tolerance = solver.number<double>(6);
    auto row = reader.take();
    while (row.words[0] != "end") {
        read_boundary(row, p);
        row = reader.take();
    }
    row.count(1);
    if (auto extra = reader.next()) extra->fail("Unexpected record after end.");
    model::validate_thermal_project(p);
    return p;
}

void write_thermal_project(std::ostream& output, const model::ThermalProject& project) {
    const auto data = canonical(project);
    output.write(data.data(), static_cast<std::streamsize>(data.size()));
    if (!output) throw std::runtime_error("Failed to write thermal project.");
}

model::ThermalProject load_thermal_project(const std::filesystem::path& path) {
    if (!std::filesystem::is_regular_file(path)) throw std::invalid_argument("Thermal project must be a regular file.");
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open thermal project.");
    return read_thermal_project(input);
}

void save_thermal_project(const std::filesystem::path& path, const model::ThermalProject& project) {
    const auto data = canonical(project);
    detail::write_atomic_file(path, [&](std::ostream& out) {
        out << data;
    }, [&](const std::filesystem::path& staged) {
        if (canonical(load_thermal_project(staged)) != data) throw std::runtime_error("Project read-back mismatch.");
    });
}

void save_thermal_result_vtk(const std::filesystem::path& path, const model::ThermalStudyResult& result) {
    detail::write_atomic_file(path, [&](std::ostream& out) {
        const auto cells = result.mesh.elements.size();
        write_legacy_vtk_conduction_3d(out, result.mesh, result.temperature.values,
            std::vector<int>(cells, 1), result.cell_fluxes, std::vector<double>(cells, result.generation_w_m3));
    });
}

} // namespace nova::io
