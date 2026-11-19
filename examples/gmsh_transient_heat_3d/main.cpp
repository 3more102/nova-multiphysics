#include "nova/fem/transient_heat_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/io/gmsh22_tetrahedral.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/io/vtk_series.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
namespace fs = std::filesystem;
using namespace nova::fem;

constexpr const char* usage =
    "Usage: nova_gmsh_transient_heat_3d mesh.msh --material VOL k [--material VOL k ...]\n"
    "       --capacity VOL C [--capacity VOL C ...] --initial T0 --dt seconds\n"
    "       --steps N [--scheme be|cn] [--material-tensor VOL xx yy zz xy xz yz]\n"
    "       [--temperature SURF T] [--flux SURF q_out] [--convection SURF h T_inf]\n"
    "       [--source Q] [--heat-source VOL Q] [--vtk final.vtk | --vtk-series result.pvd]\n"
    "       VOL and SURF accept numeric Gmsh physical tags or matching physical names.";

template<class T>
T number(std::string_view token) {
    T value{};
    const auto end = token.data() + token.size();
    const auto parsed = std::from_chars(token.data(), end, value);
    if (parsed.ec != std::errc{} || parsed.ptr != end) {
        throw std::invalid_argument("Invalid numeric argument: " + std::string(token));
    }
    if constexpr (std::is_floating_point_v<T>) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Argument must be finite: " + std::string(token));
        }
    }
    return value;
}

int physical_tag(std::string_view text, int dimension,
                 const nova::io::GmshTetrahedralMesh3D& mesh) {
    int tag{};
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), tag);
    if (parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size()) {
        if (tag < 0) throw std::invalid_argument("Physical tags must be nonnegative.");
        return tag;
    }
    for (const auto& [key, name] : mesh.physical_names) {
        if (key.first == dimension && name == text) return key.second;
    }
    throw std::invalid_argument("Unknown Gmsh physical group: " + std::string(text));
}

ScalarField3D constant(double value) {
    return [value](double, double, double) { return value; };
}

enum class BoundaryKind { temperature, flux, convection };
struct Boundary {
    BoundaryKind kind{};
    double value{};
    double ambient{};
};

struct Settings {
    std::map<int, double> materials;
    std::map<int, SymmetricTensor3D> tensor_materials;
    std::map<int, double> capacities;
    std::map<int, double> heat_sources;
    std::map<int, Boundary> boundaries;
    std::optional<double> initial;
    std::optional<double> dt;
    std::optional<std::size_t> steps;
    std::optional<double> source;
    std::optional<fs::path> vtk;
    std::optional<fs::path> vtk_series;
    TransientHeat3DScheme scheme{TransientHeat3DScheme::BackwardEuler};
    bool scheme_set{};
};

Settings read_settings(int argc, char** argv,
                       const nova::io::GmshTetrahedralMesh3D& imported) {
    Settings settings;
    int index = 2;
    const auto next = [&]() -> std::string_view {
        if (index >= argc) throw std::invalid_argument(usage);
        return argv[index++];
    };
    while (index < argc) {
        const auto option = next();
        if (option == "--material" || option == "--material-tensor") {
            const int tag = physical_tag(next(), 3, imported);
            if (option == "--material") {
                const double k = number<double>(next());
                if (!(k > 0.0) || settings.tensor_materials.contains(tag) ||
                    !settings.materials.emplace(tag, k).second) {
                    throw std::invalid_argument("Conductivity must be positive with unique volume tags.");
                }
            } else {
                const SymmetricTensor3D k{number<double>(next()), number<double>(next()),
                                          number<double>(next()), number<double>(next()),
                                          number<double>(next()), number<double>(next())};
                if (settings.materials.contains(tag) ||
                    !settings.tensor_materials.emplace(tag, k).second) {
                    throw std::invalid_argument("Duplicate tensor conductivity tag.");
                }
            }
        } else if (option == "--capacity") {
            const int tag = physical_tag(next(), 3, imported);
            const double c = number<double>(next());
            if (!(c > 0.0) || !settings.capacities.emplace(tag, c).second) {
                throw std::invalid_argument("Heat capacity must be positive with unique volume tags.");
            }
        } else if (option == "--heat-source") {
            const int tag = physical_tag(next(), 3, imported);
            if (!settings.heat_sources.emplace(tag, number<double>(next())).second) {
                throw std::invalid_argument("Duplicate volume heat source tag.");
            }
        } else if (option == "--source") {
            if (settings.source) throw std::invalid_argument("Duplicate --source.");
            settings.source = number<double>(next());
        } else if (option == "--initial") {
            if (settings.initial) throw std::invalid_argument("Duplicate --initial.");
            settings.initial = number<double>(next());
        } else if (option == "--dt") {
            if (settings.dt) throw std::invalid_argument("Duplicate --dt.");
            settings.dt = number<double>(next());
        } else if (option == "--steps") {
            if (settings.steps) throw std::invalid_argument("Duplicate --steps.");
            settings.steps = number<std::size_t>(next());
        } else if (option == "--scheme") {
            if (settings.scheme_set) throw std::invalid_argument("Duplicate --scheme.");
            settings.scheme_set = true;
            const auto scheme = next();
            if (scheme == "be") settings.scheme = TransientHeat3DScheme::BackwardEuler;
            else if (scheme == "cn") settings.scheme = TransientHeat3DScheme::CrankNicolson;
            else throw std::invalid_argument("Use --scheme be or cn.");
        } else if (option == "--vtk") {
            if (settings.vtk) throw std::invalid_argument("Duplicate --vtk.");
            const auto filename = next();
            if (filename.empty()) throw std::invalid_argument("VTK output path is empty.");
            settings.vtk = fs::path(filename);
        } else if (option == "--vtk-series") {
            if (settings.vtk_series) throw std::invalid_argument("Duplicate --vtk-series.");
            const auto filename = next();
            if (filename.empty()) throw std::invalid_argument("PVD output path is empty.");
            settings.vtk_series = fs::path(filename);
        } else if (option == "--temperature" || option == "--flux" || option == "--convection") {
            const int tag = physical_tag(next(), 2, imported);
            Boundary bc;
            bc.value = number<double>(next());
            if (option == "--temperature") bc.kind = BoundaryKind::temperature;
            else if (option == "--flux") bc.kind = BoundaryKind::flux;
            else {
                bc.kind = BoundaryKind::convection;
                if (bc.value < 0.0) throw std::invalid_argument("Convection h must be nonnegative.");
                bc.ambient = number<double>(next());
            }
            if (!settings.boundaries.emplace(tag, bc).second) {
                throw std::invalid_argument("Duplicate or conflicting boundary tag.");
            }
        } else {
            throw std::invalid_argument("Unknown argument: " + std::string(option) + "\n" + usage);
        }
    }
    if (!settings.dt || !(*settings.dt > 0.0) || !settings.steps || !(*settings.steps > 0U) ||
        !settings.initial || (settings.materials.empty() && settings.tensor_materials.empty()) ||
        settings.capacities.empty()) {
        throw std::invalid_argument(usage);
    }
    if (settings.vtk && settings.vtk_series) {
        throw std::invalid_argument("Choose either --vtk or --vtk-series, not both.");
    }
    return settings;
}

template<class T>
void require_known_tags(const std::map<int, T>& selected, const std::vector<int>& available,
                        std::string_view kind) {
    const std::set<int> known(available.begin(), available.end());
    for (const auto& [tag, value] : selected) {
        (void)value;
        if (!known.contains(tag)) {
            throw std::invalid_argument("No " + std::string(kind) +
                                        " with physical tag " + std::to_string(tag));
        }
    }
}

TransientHeat3DProblem bind_problem(nova::io::GmshTetrahedralMesh3D& imported,
                                   const Settings& settings) {
    require_known_tags(settings.materials, imported.volume_physical_tags, "volume");
    require_known_tags(settings.tensor_materials, imported.volume_physical_tags, "volume");
    require_known_tags(settings.capacities, imported.volume_physical_tags, "volume");
    require_known_tags(settings.heat_sources, imported.volume_physical_tags, "volume");
    require_known_tags(settings.boundaries, imported.exterior_physical_tags, "exterior face");

    TransientHeat3DProblem problem;
    auto& p = problem.diffusion;
    p.mesh = std::move(imported.mesh);
    if (settings.source) p.source = constant(*settings.source);
    const bool tensor_mode = !settings.tensor_materials.empty();
    const std::size_t cells = imported.volume_physical_tags.size();
    problem.cell_volumetric_capacity.reserve(cells);
    p.cell_source.reserve(cells);
    if (tensor_mode) p.cell_diffusivity_tensors.reserve(cells);
    else p.cell_diffusivity.reserve(cells);
    for (const int tag : imported.volume_physical_tags) {
        const auto scalar = settings.materials.find(tag);
        const auto tensor = settings.tensor_materials.find(tag);
        const auto capacity = settings.capacities.find(tag);
        if ((scalar == settings.materials.end() && tensor == settings.tensor_materials.end()) ||
            capacity == settings.capacities.end()) {
            throw std::invalid_argument("Missing conductivity or capacity for volume tag " +
                                        std::to_string(tag));
        }
        if (tensor_mode) {
            p.cell_diffusivity_tensors.push_back(tensor != settings.tensor_materials.end()
                ? tensor->second : SymmetricTensor3D{scalar->second, scalar->second, scalar->second});
        } else {
            p.cell_diffusivity.push_back(scalar->second);
        }
        problem.cell_volumetric_capacity.push_back(capacity->second);
        const auto source = settings.heat_sources.find(tag);
        p.cell_source.push_back(source == settings.heat_sources.end() ? 0.0 : source->second);
    }
    for (std::size_t i = 0; i < imported.exterior_faces.size(); ++i) {
        const auto found = settings.boundaries.find(imported.exterior_physical_tags[i]);
        if (found == settings.boundaries.end()) continue;
        const auto& face = imported.exterior_faces[i];
        const auto& bc = found->second;
        switch (bc.kind) {
        case BoundaryKind::temperature: p.dirichlet_faces.push_back({face, constant(bc.value)}); break;
        case BoundaryKind::flux: p.neumann_faces.push_back({face, constant(bc.value)}); break;
        case BoundaryKind::convection:
            p.robin_faces.push_back({face, constant(bc.value), constant(bc.ambient)}); break;
        }
    }
    problem.initial_temperature.assign(p.mesh.nodes.size(), *settings.initial);
    return problem;
}

void export_vtk(const fs::path& target, const TransientHeat3DProblem& problem,
                const std::vector<double>& temperature, const std::vector<int>& volume_tags,
                const Settings& settings) {
    const auto& p = problem.diffusion;
    const auto flux = evaluate_p1_tetrahedral_flux_3d(p, temperature);
    std::vector<double> generation;
    generation.reserve(volume_tags.size());
    for (const int tag : volume_tags) {
        const auto it = settings.heat_sources.find(tag);
        const double value = settings.source.value_or(0.0) +
            (it == settings.heat_sources.end() ? 0.0 : it->second);
        if (!std::isfinite(value)) throw std::overflow_error("Effective source overflow.");
        generation.push_back(value);
    }
    fs::path temporary;
    bool created = false;
    for (std::size_t attempt = 0; attempt < 1000U; ++attempt) {
        temporary = target;
        temporary += ".tmp." + std::to_string(attempt);
        if (fs::create_directory(temporary)) {
            created = true;
            break;
        }
    }
    if (!created) throw std::runtime_error("Cannot create a temporary VTK result directory.");
    const auto frame = temporary / "frame.vtk";
    struct Cleanup {
        const fs::path& frame;
        const fs::path& directory;
        ~Cleanup() {
            std::error_code ignored;
            fs::remove(frame, ignored);
            fs::remove(directory, ignored);
        }
    } cleanup{frame, temporary};
    std::ofstream output(frame);
    if (!output) throw std::runtime_error("Cannot open VTK result output.");
    nova::io::write_legacy_vtk_conduction_3d(
        output, p.mesh, temperature, volume_tags, flux, generation);
    output.close();
    if (!output) throw std::runtime_error("Unable to complete VTK result output.");
    fs::rename(frame, target);
}


std::vector<double> make_heat_generation(
    const std::vector<int>& volume_tags, const Settings& settings
) {
    std::vector<double> generation;
    generation.reserve(volume_tags.size());
    for (const int tag : volume_tags) {
        const auto found = settings.heat_sources.find(tag);
        const double value = settings.source.value_or(0.0) +
            (found == settings.heat_sources.end() ? 0.0 : found->second);
        if (!std::isfinite(value)) {
            throw std::overflow_error("Effective cell heat generation overflowed.");
        }
        generation.push_back(value);
    }
    return generation;
}

void write_series_frame(
    const fs::path& target, const TransientHeat3DProblem& problem,
    const std::vector<double>& temperature, const std::vector<int>& volume_tags,
    const std::vector<double>& generation
) {
    const auto flux = evaluate_p1_tetrahedral_flux_3d(problem.diffusion, temperature);
    fs::path temporary;
    bool created = false;
    for (std::size_t attempt = 0; attempt < 1000U; ++attempt) {
        temporary = target;
        temporary += ".tmp." + std::to_string(attempt);
        if (fs::create_directory(temporary)) {
            created = true;
            break;
        }
    }
    if (!created) throw std::runtime_error("Cannot create a temporary VTK frame directory.");
    const auto staged = temporary / "frame.vtk";
    struct Cleanup {
        const fs::path& file;
        const fs::path& directory;
        ~Cleanup() {
            std::error_code ignored;
            fs::remove(file, ignored);
            fs::remove(directory, ignored);
        }
    } cleanup{staged, temporary};
    std::ofstream output(staged);
    if (!output) throw std::runtime_error("Cannot open staged VTK frame.");
    nova::io::write_legacy_vtk_conduction_3d(
        output, problem.diffusion.mesh, temperature, volume_tags, flux, generation);
    output.close();
    if (!output) throw std::runtime_error("Cannot finish staged VTK frame.");
    fs::rename(staged, target);
}

void write_series_collection(
    const fs::path& destination,
    const std::vector<nova::io::VtkTimeSeriesEntry>& frames
) {
    fs::path temporary;
    bool created = false;
    for (std::size_t attempt = 0; attempt < 1000U; ++attempt) {
        temporary = destination;
        temporary += ".tmp." + std::to_string(attempt);
        if (fs::create_directory(temporary)) {
            created = true;
            break;
        }
    }
    if (!created) throw std::runtime_error("Cannot create a temporary PVD directory.");
    const auto staged = temporary / "collection.pvd";
    struct Cleanup {
        const fs::path& file;
        const fs::path& directory;
        ~Cleanup() {
            std::error_code ignored;
            fs::remove(file, ignored);
            fs::remove(directory, ignored);
        }
    } cleanup{staged, temporary};
    std::ofstream output(staged);
    if (!output) throw std::runtime_error("Cannot open staged PVD collection.");
    nova::io::write_pvd_collection(output, frames);
    output.close();
    if (!output) throw std::runtime_error("Cannot finish staged PVD collection.");
    fs::rename(staged, destination);
}

std::string series_prefix(const fs::path& destination) {
    if (destination.extension() != ".pvd") {
        throw std::invalid_argument("--vtk-series requires a .pvd destination.");
    }
    const std::string prefix = destination.stem().string();
    if (prefix.empty() || prefix.front() == '.' ||
        prefix.find("..") != std::string::npos ||
        !std::all_of(prefix.begin(), prefix.end(), [](char c) {
            return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                   (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.';
        })) {
        throw std::invalid_argument("Use a safe ASCII .pvd destination basename.");
    }
    return prefix;
}

std::string series_frame_name(const std::string& prefix, std::size_t step) {
    std::ostringstream name;
    name << prefix << '_' << std::setw(6) << std::setfill('0') << step << ".vtk";
    return name.str();
}

void validate_series_destination(
    const fs::path& destination, const std::string& prefix, std::size_t steps
) {
    if (steps == std::numeric_limits<std::size_t>::max()) {
        throw std::length_error("Transient step count is too large for a PVD series.");
    }
    const fs::path parent = destination.parent_path().empty()
        ? fs::path(".") : destination.parent_path();
    if (!fs::is_directory(parent)) {
        throw std::invalid_argument("PVD output directory does not exist.");
    }
    if (fs::exists(destination)) {
        throw std::invalid_argument("Refusing to overwrite an existing PVD collection.");
    }
    for (std::size_t step = 0; step <= steps; ++step) {
        const auto frame = parent / series_frame_name(prefix, step);
        if (fs::exists(frame)) {
            throw std::invalid_argument("Refusing to overwrite existing VTK frame: " +
                                        frame.string());
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3) throw std::invalid_argument(usage);
        std::ifstream file(argv[1]);
        if (!file) throw std::runtime_error("Cannot open Gmsh mesh: " + std::string(argv[1]));
        auto imported = nova::io::read_gmsh22_tetrahedra(file);
        const auto settings = read_settings(argc, argv, imported);
        if (settings.vtk) {
            const fs::path parent = settings.vtk->parent_path().empty()
                ? fs::path(".") : settings.vtk->parent_path();
            if (!fs::is_directory(parent)) {
                throw std::invalid_argument("VTK output directory does not exist.");
            }
            if (fs::exists(*settings.vtk)) {
                if (fs::equivalent(*settings.vtk, fs::path(argv[1]))) {
                    throw std::invalid_argument("VTK path must not overwrite the input Gmsh mesh.");
                }
                throw std::invalid_argument("Refusing to overwrite an existing VTK result.");
            }
        }
        std::string prefix;
        if (settings.vtk_series) {
            prefix = series_prefix(*settings.vtk_series);
            validate_series_destination(*settings.vtk_series, prefix, *settings.steps);
        }
        auto problem = bind_problem(imported, settings);
        const auto generation = make_heat_generation(imported.volume_physical_tags, settings);
        TransientHeat3DOptions options;
        options.scheme = settings.scheme;
        options.time_step = *settings.dt;
        options.steps = *settings.steps;
        options.solver_options.relative_tolerance = 1e-11;
        options.solver_options.absolute_tolerance = 0.0;
        options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        options.retain_step_solutions = false; // The CLI only needs the final field and diagnostics.
        std::vector<nova::io::VtkTimeSeriesEntry> frames;
        if (settings.vtk_series) {
            frames.reserve(*settings.steps + 1U);
            options.on_step = [&](std::size_t step, double time,
                                  const std::vector<double>& temperature) {
                const auto name = series_frame_name(prefix, step);
                const auto frame = settings.vtk_series->parent_path() / name;
                write_series_frame(frame, problem, temperature,
                                   imported.volume_physical_tags, generation);
                frames.push_back({time, name});
            };
        }
        const auto result = solve_transient_heat_3d(problem, options);
        if (settings.vtk) {
            export_vtk(*settings.vtk, problem, result.temperature,
                       imported.volume_physical_tags, settings);
        }
        if (settings.vtk_series) {
            write_series_collection(*settings.vtk_series, frames);
        }
        const auto [min_temp, max_temp] = std::minmax_element(
            result.temperature.begin(), result.temperature.end());
        double max_balance = 0.0;
        for (const auto& step : result.diagnostics) {
            max_balance = std::max(max_balance, std::abs(step.balance_residual));
        }
        std::cout << std::setprecision(17)
                  << "NOVA Gmsh transient 3D heat study\n"
                  << "steps " << options.steps << " final_time " << result.final_time
                  << " min_temperature " << *min_temp << " max_temperature " << *max_temp << "\n"
                  << "max_abs_energy_balance_residual " << max_balance << "\n";
        if (settings.vtk_series) {
            std::cout << "ParaView collection " << settings.vtk_series->string()
                      << " frames " << frames.size() << "\n";
        }
        if (settings.vtk) {
            std::cout << "VTK result " << settings.vtk->string() << "\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
