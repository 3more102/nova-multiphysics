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
#include <sstream>
#include <system_error>
#include <map>
#include <optional>
#include <set>
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
    "       [--start-time seconds] [--temperature SURF T] [--flux SURF q_out]\n"
    "       [--convection SURF h T_inf] [--source Q] [--heat-source VOL Q]\n"
    "       [--vtk final.vtk | --vtk-series result.pvd]\n"
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
    std::optional<double> start_time;
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
        } else if (option == "--start-time") {
            if (settings.start_time) throw std::invalid_argument("Duplicate --start-time.");
            settings.start_time = number<double>(next());
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
    p.source = constant(settings.source.value_or(0.0));
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

std::vector<double> make_heat_generation(
    const std::vector<int>& volume_tags, const Settings& settings
) {
    std::vector<double> generation;
    generation.reserve(volume_tags.size());
    for (const int tag : volume_tags) {
        const auto it = settings.heat_sources.find(tag);
        const double value = settings.source.value_or(0.0) +
            (it == settings.heat_sources.end() ? 0.0 : it->second);
        if (!std::isfinite(value)) throw std::overflow_error("Effective source overflow.");
        generation.push_back(value);
    }
    return generation;
}

void write_vtk_frame(
    const fs::path& target, const TransientHeat3DProblem& problem,
    const std::vector<double>& temperature, const std::vector<int>& volume_tags,
    const std::vector<double>& generation
) {
    const auto flux = evaluate_p1_tetrahedral_flux_3d(problem.diffusion, temperature);
    std::ofstream output(target);
    if (!output) throw std::runtime_error("Cannot open VTK frame: " + target.string());
    nova::io::write_legacy_vtk_conduction_3d(
        output, problem.diffusion.mesh, temperature, volume_tags, flux, generation);
    output.close();
    if (!output) throw std::runtime_error("Cannot finish VTK frame: " + target.string());
}

void export_vtk(
    const fs::path& target, const TransientHeat3DProblem& problem,
    const std::vector<double>& temperature, const std::vector<int>& volume_tags,
    const std::vector<double>& generation
) {
    if (fs::exists(target)) {
        throw std::invalid_argument("VTK output already exists: " + target.string());
    }
    fs::path temporary;
    for (std::size_t attempt = 0; ; ++attempt) {
        temporary = target;
        temporary += ".tmp." + std::to_string(attempt);
        if (fs::create_directory(temporary)) break;
    }
    struct Cleanup {
        const fs::path& directory;
        ~Cleanup() {
            std::error_code ignored;
            fs::remove_all(directory, ignored);
        }
    } cleanup{temporary};
    const auto frame = temporary / "frame.vtk";
    write_vtk_frame(frame, problem, temperature, volume_tags, generation);
    fs::rename(frame, target);
}

class VtkSeriesWriter {
public:
    VtkSeriesWriter(
        fs::path destination, const TransientHeat3DProblem& problem,
        const std::vector<int>& volume_tags, const std::vector<double>& generation,
        std::size_t expected_frames
    )
        : destination_(std::move(destination)),
          parent_(destination_.parent_path().empty() ? fs::path(".") : destination_.parent_path()),
          prefix_(destination_.stem().string()),
          problem_(problem),
          volume_tags_(volume_tags),
          generation_(generation),
          expected_frames_(expected_frames) {
        if (destination_.extension() != ".pvd" || !valid_prefix(prefix_)) {
            throw std::invalid_argument("Use a .pvd path with a safe ASCII basename.");
        }
        if (!fs::is_directory(parent_)) {
            throw std::invalid_argument("PVD output directory does not exist: " + parent_.string());
        }
        if (expected_frames_ == 0U) {
            throw std::invalid_argument("PVD series requires at least one frame.");
        }
        if (fs::exists(destination_)) {
            throw std::invalid_argument("PVD output already exists: " + destination_.string());
        }
        // Fail before the expensive solve if any deterministic frame name collides.
        for (std::size_t step = 0; step < expected_frames_; ++step) {
            const auto frame = parent_ / frame_filename(step);
            if (fs::exists(frame)) {
                throw std::invalid_argument("VTK frame already exists: " + frame.string());
            }
        }
        for (std::size_t attempt = 0; ; ++attempt) {
            staging_directory_ = parent_ /
                ("." + prefix_ + ".nova-tmp-" + std::to_string(attempt));
            if (fs::create_directory(staging_directory_)) break;
            if (attempt == std::numeric_limits<std::size_t>::max()) {
                throw std::runtime_error("Unable to allocate a temporary PVD output directory.");
            }
        }
    }

    VtkSeriesWriter(const VtkSeriesWriter&) = delete;
    VtkSeriesWriter& operator=(const VtkSeriesWriter&) = delete;

    ~VtkSeriesWriter() noexcept {
        if (!staging_directory_.empty()) {
            std::error_code ignored;
            fs::remove_all(staging_directory_, ignored);
        }
    }

    void write_step(
        std::size_t step, double time, const std::vector<double>& temperature
    ) {
        if (step != entries_.size() || step >= expected_frames_) {
            throw std::logic_error("PVD observer steps must be consecutive and start at zero.");
        }
        if (!std::isfinite(time)) {
            throw std::invalid_argument("PVD frame time must be finite.");
        }
        const auto filename = frame_filename(step);
        write_vtk_frame(staging_directory_ / filename, problem_, temperature,
                        volume_tags_, generation_);
        entries_.push_back({time, filename});
    }

    void commit() {
        if (entries_.size() != expected_frames_) {
            throw std::logic_error("Cannot publish an incomplete PVD time series.");
        }
        const auto staged_collection = staging_directory_ / destination_.filename();
        {
            std::ofstream output(staged_collection);
            if (!output) throw std::runtime_error("Cannot open staged PVD collection.");
            nova::io::write_pvd_collection(output, entries_);
            output.close();
            if (!output) throw std::runtime_error("Cannot finish staged PVD collection.");
        }

        std::vector<fs::path> published;
        published.reserve(entries_.size() + 1U);
        try {
            for (const auto& entry : entries_) {
                const auto final_frame = parent_ / entry.filename;
                if (fs::exists(final_frame)) {
                    throw std::runtime_error("VTK frame appeared during export: " + final_frame.string());
                }
                fs::rename(staging_directory_ / entry.filename, final_frame);
                published.push_back(final_frame);
            }
            if (fs::exists(destination_)) {
                throw std::runtime_error("PVD collection appeared during export: " + destination_.string());
            }
            // Publish the manifest last: readers never see a collection before its frames.
            fs::rename(staged_collection, destination_);
            published.push_back(destination_);
        } catch (...) {
            for (const auto& path : published) {
                std::error_code ignored;
                fs::remove(path, ignored);
            }
            throw;
        }

        std::error_code ignored;
        fs::remove_all(staging_directory_, ignored);
        if (!ignored) staging_directory_.clear();
    }

    [[nodiscard]] std::size_t frames_written() const noexcept {
        return entries_.size();
    }

private:
    [[nodiscard]] static bool valid_prefix(const std::string& prefix) noexcept {
        if (prefix.empty() || prefix.front() == '.' ||
            prefix.find("..") != std::string::npos) {
            return false;
        }
        for (const char c : prefix) {
            const bool alpha = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
            const bool digit = c >= '0' && c <= '9';
            if (!alpha && !digit && c != '_' && c != '-' && c != '.') return false;
        }
        return true;
    }

    [[nodiscard]] std::string frame_filename(std::size_t step) const {
        std::ostringstream name;
        name << prefix_ << '_' << std::setw(6) << std::setfill('0') << step << ".vtk";
        return name.str();
    }

    fs::path destination_;
    fs::path parent_;
    std::string prefix_;
    const TransientHeat3DProblem& problem_;
    const std::vector<int>& volume_tags_;
    const std::vector<double>& generation_;
    std::size_t expected_frames_{};
    fs::path staging_directory_;
    std::vector<nova::io::VtkTimeSeriesEntry> entries_;
};
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3) throw std::invalid_argument(usage);
        const fs::path mesh_path(argv[1]);
        std::ifstream file(mesh_path);
        if (!file) throw std::runtime_error("Cannot open Gmsh mesh: " + mesh_path.string());
        auto imported = nova::io::read_gmsh22_tetrahedra(file);
        const auto settings = read_settings(argc, argv, imported);

        const auto reject_existing_output = [&](const std::optional<fs::path>& path) {
            if (!path || !fs::exists(*path)) return;
            if (fs::equivalent(*path, mesh_path)) {
                throw std::invalid_argument("Output path must not overwrite the input Gmsh mesh.");
            }
            throw std::invalid_argument("Output path already exists: " + path->string());
        };
        reject_existing_output(settings.vtk);
        reject_existing_output(settings.vtk_series);
        if (settings.vtk && settings.vtk->extension() != ".vtk") {
            throw std::invalid_argument("--vtk output must use the .vtk extension.");
        }
        if (settings.vtk) {
            const auto parent = settings.vtk->parent_path().empty()
                ? fs::path(".") : settings.vtk->parent_path();
            if (!fs::is_directory(parent)) {
                throw std::invalid_argument("VTK output directory does not exist: " + parent.string());
            }
        }

        auto problem = bind_problem(imported, settings);
        const auto generation = make_heat_generation(imported.volume_physical_tags, settings);
        TransientHeat3DOptions options;
        options.scheme = settings.scheme;
        options.time_step = *settings.dt;
        options.start_time = settings.start_time.value_or(0.0);
        options.steps = *settings.steps;
        options.retain_step_solutions = false;
        options.solver_options.relative_tolerance = 1e-11;
        options.solver_options.absolute_tolerance = 0.0;
        options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;

        std::optional<VtkSeriesWriter> series;
        if (settings.vtk_series) {
            if (options.steps == std::numeric_limits<std::size_t>::max()) {
                throw std::length_error("Too many steps for a PVD series.");
            }
            series.emplace(*settings.vtk_series, problem, imported.volume_physical_tags,
                           generation, options.steps + 1U);
            options.on_step = [&series](
                std::size_t step, double time, const std::vector<double>& temperature
            ) {
                if (!series) throw std::logic_error("PVD observer has no output writer.");
                series->write_step(step, time, temperature);
            };
        }

        const auto result = solve_transient_heat_3d(problem, options);
        if (series) {
            series->commit();
        } else if (settings.vtk) {
            export_vtk(*settings.vtk, problem, result.temperature,
                       imported.volume_physical_tags, generation);
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
        if (series) {
            std::cout << "ParaView collection: " << settings.vtk_series->string()
                      << " (" << series->frames_written() << " frames)\n";
        } else if (settings.vtk) {
            std::cout << "VTK result: " << settings.vtk->string() << "\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
