#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/energy_balance_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"
#include "nova/fem/tagged_boundary_flux_3d.hpp"
#include "nova/io/gmsh_tetrahedral.hpp"
#include "nova/io/vtk_legacy_3d.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

namespace fs = std::filesystem;
using namespace nova::fem;

constexpr const char* usage =
    "Usage: nova_gmsh_diffusion_3d mesh.msh <dirichlet-tag-or-name>\n"
    "   or: nova_gmsh_diffusion_3d mesh.msh --material volume-tag-or-name k ...\n"
    "       [--temperature face-tag-or-name T] [--flux face-tag-or-name q_out]\n"
    "       [--material-tensor volume-tag-or-name xx yy zz xy xz yz]\n"
    "       [--convection face-tag-or-name h T_ambient] [--source Q]\n"
    "       [--heat-source volume-tag-or-name Q] ... [--vtk result.vtk]\n"
    "       [--flux-report] [--balance-report]";

template<class T>
T number(std::string_view token) {
    T value{};
    const auto end = token.data() + token.size();
    const auto parsed = std::from_chars(token.data(), end, value);
    if (parsed.ec != std::errc{} || parsed.ptr != end || !std::isfinite(value)) {
        throw std::invalid_argument("Invalid finite number: " + std::string(token));
    }
    return value;
}

int physical_tag(std::string_view token, int dimension,
                 const nova::io::GmshTetrahedralMesh3D& mesh) {
    int numeric{};
    const auto parsed = std::from_chars(token.data(), token.data() + token.size(), numeric);
    if (parsed.ec == std::errc{} && parsed.ptr == token.data() + token.size()) {
        if (numeric < 0) throw std::invalid_argument("Physical tags must be nonnegative.");
        return numeric;
    }
    for (const auto& [key, name] : mesh.physical_names) {
        if (key.first == dimension && name == token) return key.second;
    }
    throw std::invalid_argument("Unknown Gmsh physical group: " + std::string(token));
}

ScalarField3D constant(double value) {
    return [value](double, double, double) { return value; };
}

enum class BoundaryKind { temperature, flux, convection };
struct Boundary {
    BoundaryKind kind;
    double value;
    double ambient{};
};

struct Settings {
    std::map<int, double> materials;
    std::map<int, SymmetricTensor3D> tensor_materials;
    std::map<int, double> heat_sources;
    std::map<int, Boundary> surfaces;
    std::optional<double> source;
    std::optional<fs::path> vtk;
    bool unit_example{};
    bool flux_report{};
    bool balance_report{};
};

Settings read_settings(int argc, char** argv,
                       const nova::io::GmshTetrahedralMesh3D& mesh) {
    if (argc < 3) throw std::invalid_argument(usage);
    Settings settings;
    if (argc == 3 && !std::string_view(argv[2]).starts_with("--")) {
        const int tag = physical_tag(argv[2], 2, mesh);
        if (tag == 0) throw std::invalid_argument("The example needs a positive Dirichlet tag.");
        settings.unit_example = true;
        settings.source = 1.0;
        settings.surfaces.emplace(tag, Boundary{BoundaryKind::temperature, 0.0});
        return settings;
    }

    int index = 2;
    auto next = [&]() -> std::string_view {
        if (index == argc) throw std::invalid_argument(usage);
        return argv[index++];
    };
    while (index < argc) {
        const auto option = next();
        if (option == "--material") {
            const int tag = physical_tag(next(), 3, mesh);
            const double k = number<double>(next());
            if (!(k > 0.0)) throw std::invalid_argument("Material conductivity must be positive.");
            if (settings.tensor_materials.contains(tag) || !settings.materials.emplace(tag, k).second) {
                throw std::invalid_argument("Duplicate material tag " + std::to_string(tag));
            }
        } else if (option == "--material-tensor") {
            const int tag = physical_tag(next(), 3, mesh);
            // Components: xx yy zz xy xz yz, global Cartesian axes.
            const SymmetricTensor3D tensor{number<double>(next()), number<double>(next()),
                                           number<double>(next()), number<double>(next()),
                                           number<double>(next()), number<double>(next())};
            if (settings.materials.contains(tag) ||
                !settings.tensor_materials.emplace(tag, tensor).second) {
                throw std::invalid_argument("Duplicate material tag " + std::to_string(tag));
            }
        } else if (option == "--heat-source") {
            const int tag = physical_tag(next(), 3, mesh);
            const double heat = number<double>(next());
            if (!settings.heat_sources.emplace(tag, heat).second) {
                throw std::invalid_argument("Duplicate heat source for volume tag " + std::to_string(tag));
            }
        } else if (option == "--source") {
            if (settings.source) throw std::invalid_argument("Duplicate --source.");
            settings.source = number<double>(next());
        } else if (option == "--vtk") {
            if (settings.vtk) throw std::invalid_argument("Duplicate --vtk.");
            settings.vtk = fs::path(next());
        } else if (option == "--flux-report") {
            if (settings.flux_report) throw std::invalid_argument("Duplicate --flux-report.");
            settings.flux_report = true;
        } else if (option == "--balance-report") {
            if (settings.balance_report) throw std::invalid_argument("Duplicate --balance-report.");
            settings.balance_report = true;
        } else if (option == "--temperature" || option == "--flux" || option == "--convection") {
            const int tag = physical_tag(next(), 2, mesh);
            Boundary condition{BoundaryKind::temperature, number<double>(next())};
            if (option == "--flux") condition.kind = BoundaryKind::flux;
            if (option == "--convection") {
                condition.kind = BoundaryKind::convection;
                if (condition.value < 0.0) throw std::invalid_argument("Convection h must be nonnegative.");
                condition.ambient = number<double>(next());
            }
            if (!settings.surfaces.emplace(tag, condition).second) {
                throw std::invalid_argument("Duplicate surface condition on tag " + std::to_string(tag));
            }
        } else {
            throw std::invalid_argument("Unknown option " + std::string(option) + "\n" + usage);
        }
    }
    if (settings.materials.empty() && settings.tensor_materials.empty()) {
        throw std::invalid_argument("Supply at least one --material or --material-tensor.");
    }
    return settings;
}

template<class Value>
void check_selected_tags(const std::map<int, Value>& selected, const std::vector<int>& mesh_tags,
                         std::string_view entity) {
    const std::set<int> present(mesh_tags.begin(), mesh_tags.end());
    for (const auto& entry : selected) {
        if (!present.contains(entry.first)) {
            throw std::invalid_argument("No " + std::string(entity) + " with physical tag " +
                                        std::to_string(entry.first));
        }
    }
}

ScalarDiffusion3DProblem bind_problem(nova::io::GmshTetrahedralMesh3D& imported,
                                      const Settings& settings) {
    ScalarDiffusion3DProblem problem;
    problem.mesh = std::move(imported.mesh);
    problem.source = constant(settings.source.value_or(0.0));
    if (settings.unit_example) {
        problem.diffusivity = constant(1.0);
    } else {
        check_selected_tags(settings.materials, imported.volume_physical_tags, "volume");
        check_selected_tags(settings.tensor_materials, imported.volume_physical_tags, "volume");
        const bool any_tensor = !settings.tensor_materials.empty();
        if (any_tensor) {
            problem.cell_diffusivity_tensors.reserve(imported.volume_physical_tags.size());
        } else {
            problem.cell_diffusivity.reserve(imported.volume_physical_tags.size());
        }
        for (int tag : imported.volume_physical_tags) {
            const auto scalar = settings.materials.find(tag);
            const auto tensor = settings.tensor_materials.find(tag);
            if (scalar == settings.materials.end() &&
                tensor == settings.tensor_materials.end()) {
                throw std::invalid_argument("Missing conductivity for volume tag " + std::to_string(tag));
            }
            if (any_tensor) {
                // Promote isotropic regions to diagonal tensors for mixed media.
                const SymmetricTensor3D k = tensor == settings.tensor_materials.end()
                    ? SymmetricTensor3D{scalar->second, scalar->second, scalar->second}
                    : tensor->second;
                problem.cell_diffusivity_tensors.push_back(k);
            } else {
                problem.cell_diffusivity.push_back(scalar->second);
            }
        }
    }
    check_selected_tags(settings.heat_sources, imported.volume_physical_tags, "volume");
    if (!settings.heat_sources.empty()) {
        problem.cell_source.reserve(imported.volume_physical_tags.size());
        for (const int tag : imported.volume_physical_tags) {
            const auto found = settings.heat_sources.find(tag);
            problem.cell_source.push_back(
                found == settings.heat_sources.end() ? 0.0 : found->second);
        }
    }
    check_selected_tags(settings.surfaces, imported.exterior_physical_tags, "exterior face");
    for (std::size_t i = 0; i < imported.exterior_faces.size(); ++i) {
        const auto found = settings.surfaces.find(imported.exterior_physical_tags[i]);
        if (found == settings.surfaces.end()) continue;
        const auto& [kind, value, ambient] = found->second;
        const auto& face = imported.exterior_faces[i];
        switch (kind) {
        case BoundaryKind::temperature: problem.dirichlet_faces.push_back({face, constant(value)}); break;
        case BoundaryKind::flux: problem.neumann_faces.push_back({face, constant(value)}); break;
        case BoundaryKind::convection:
            problem.robin_faces.push_back({face, constant(value), constant(ambient)}); break;
        }
    }
    return problem;
}

// Publish only a completely written frame. An exclusively created sibling
// directory also separates simultaneous writers and keeps rename on one volume.
void export_result(const fs::path& target, const ScalarDiffusion3DProblem& problem,
                   const std::vector<double>& temperature, const std::vector<int>& tags,
                   double uniform_source) {
    const auto fields = evaluate_p1_tetrahedral_flux_3d(problem, temperature);
    std::vector<double> effective_generation = problem.cell_source;
    for (double& value : effective_generation) {
        value += uniform_source;
        if (!std::isfinite(value)) {
            throw std::overflow_error("Effective cell heat generation is non-finite.");
        }
    }
    fs::path temporary;
    for (std::size_t attempt = 0; ; ++attempt) {
        temporary = target;
        temporary += ".tmp." + std::to_string(attempt);
        if (fs::create_directory(temporary)) break;
    }
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
    if (!output) throw std::runtime_error("Cannot open VTK output: " + frame.string());
    nova::io::write_legacy_vtk_conduction_3d(
        output, problem.mesh, temperature, tags, fields, effective_generation);
    output.close();
    if (!output) throw std::runtime_error("Could not finish writing VTK output.");
    fs::rename(frame, target);
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3) throw std::invalid_argument(usage);
        std::ifstream input(argv[1]);
        if (!input) throw std::runtime_error("Cannot open Gmsh input: " + std::string(argv[1]));
        auto imported = nova::io::read_gmsh_tetrahedra(input);
        const auto settings = read_settings(argc, argv, imported);
        if (settings.vtk && fs::exists(*settings.vtk) && fs::equivalent(*settings.vtk, argv[1])) {
            throw std::invalid_argument("VTK output would overwrite the input mesh.");
        }
        auto problem = bind_problem(imported, settings);

        nova::linalg::SolverOptions options;
        options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        options.relative_tolerance = settings.unit_example ? 1e-11 : 1e-12;
        if (!settings.unit_example) options.absolute_tolerance = 0.0;
        const auto result = solve_scalar_diffusion_3d(problem, options);
        if (!result.linear_solver.converged) throw std::runtime_error(result.linear_solver.message);
        if (settings.vtk) export_result(
            *settings.vtk, problem, result.values, imported.volume_physical_tags,
            settings.source.value_or(0.0));
        if (settings.flux_report) {
            const auto fields = evaluate_p1_tetrahedral_flux_3d(problem, result.values);
            const auto report = evaluate_tagged_boundary_flux_3d(
                problem.mesh, imported.exterior_faces, imported.exterior_physical_tags, fields);
            long double net_rate = 0;
            std::cout << std::setprecision(17)
                      << "Boundary conductive flux (raw P1, non-equilibrated):\n";
            for (const auto& group : report) {
                std::cout << "  physical tag " << group.physical_tag
                          << ": faces=" << group.face_count
                          << " area_m2=" << group.area
                          << " outward_W=" << group.outward_heat_rate << '\n';
                net_rate += group.outward_heat_rate;
            }
            if (!std::isfinite(net_rate) ||
                std::abs(net_rate) > static_cast<long double>(std::numeric_limits<double>::max())) {
                throw std::overflow_error("Net outward heat rate is not finite.");
            }
            std::cout << "Net outward conductive heat rate [W]: "
                      << static_cast<double>(net_rate) << '\n';
        }
        if (settings.balance_report) {
            const auto balance = evaluate_diffusion_energy_balance_3d(problem, result.values);
            std::cout << std::setprecision(17)
                      << "FEM generated heat [W]: " << balance.generated_heat << '\n'
                      << "FEM Neumann outward heat [W]: " << balance.neumann_outward_heat << '\n'
                      << "FEM Robin outward heat [W]: " << balance.robin_outward_heat << '\n'
                      << "FEM constrained-node outward reaction [W]: "
                      << balance.constraint_outward_heat << '\n'
                      << "FEM signed energy balance defect [W]: "
                      << balance.signed_balance_defect << '\n'
                      << "FEM maximum free-equation residual [W]: "
                      << balance.max_free_equation_residual << '\n'
                      << "FEM constrained nodes: " << balance.constrained_nodes << '\n'
                      << "FEM free nodes: " << balance.free_nodes << '\n';
        }
        const auto [lo, hi] = std::minmax_element(result.values.begin(), result.values.end());
        std::cout << std::setprecision(17)
                  << "NOVA imported Gmsh P1 tetrahedral diffusion\n"
                  << "Nodes: " << problem.mesh.nodes.size() << '\n'
                  << "Tetrahedra: " << problem.mesh.elements.size() << '\n'
                  << "Dirichlet faces: " << problem.dirichlet_faces.size() << '\n'
                  << "CG iterations: " << result.linear_solver.iterations << '\n'
                  << "True residual: " << result.linear_solver.residual_norm << '\n'
                  << "Scalar range: " << *lo << " to " << *hi << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
