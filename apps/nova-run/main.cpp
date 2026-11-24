#include "nova/io/thermal_project.hpp"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <span>
#include <vector>
#include <string_view>

namespace {
namespace fs = std::filesystem;

constexpr std::string_view usage =
    "Usage: nova-run check PROJECT\n"
    "       nova-run show PROJECT\n"
    "       nova-run save-as PROJECT OUTPUT.nova\n"
    "       nova-run solve PROJECT [--vtk OUTPUT.vtk]\n";

int run(std::span<const fs::path> args) {
    const auto argc = args.size();
    if (argc == 2 && args[1] == "--help") {
        std::cout << usage;
        return 0;
    }
    if (argc < 3) {
        std::cerr << usage;
        return 2;
    }
    const auto& command = args[1];
    const bool check = command == "check" && argc == 3;
    const bool show = command == "show" && argc == 3;
    const bool save = command == "save-as" && argc == 4;
    const bool export_vtk = command == "solve" && argc == 5 && args[3] == "--vtk";
    if (!check && !show && !save && !export_vtk && !(command == "solve" && argc == 3)) {
        std::cerr << usage;
        return 2;
    }
    const auto& input = args[2];
    const auto project = nova::io::load_thermal_project(input);
    if (show) {
        nova::io::write_thermal_project(std::cout, project);
        return 0;
    }
    if (save) {
        const auto& output = args[3];
        if (output.extension() != ".nova") throw std::invalid_argument("Save destination requires .nova.");
        nova::io::save_thermal_project(output, project);
        return 0;
    }
    if (check) {
        const auto problem = nova::model::prepare_thermal_project(project);
        std::cout << "NOVA_CHECK 1\nproject " << project.id << "\nnodes " << problem.mesh.nodes.size()
                  << "\ntetrahedra " << problem.mesh.elements.size() << '\n';
        return 0;
    }
    fs::path output;
    if (export_vtk) {
        output = args[4];
        if (output.extension() != ".vtk" || (fs::exists(output) && fs::equivalent(input, output))) {
            throw std::invalid_argument("VTK destination must be a .vtk file distinct from the input project.");
        }
    }
    const auto result = nova::model::solve_thermal_project(project);
    if (export_vtk) nova::io::save_thermal_result_vtk(output, result);
    const auto [cold, hot] = std::minmax_element(result.temperature.values.begin(), result.temperature.values.end());
    const auto& cg = result.temperature.linear_solver;
    const auto& balance = result.balance;
    std::cout.imbue(std::locale::classic());
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "NOVA_RESULT 1\nproject " << project.id << "\nstudy " << project.study_id
              << "\nnodes " << result.mesh.nodes.size() << "\ntetrahedra " << result.mesh.elements.size()
              << "\niterations " << cg.iterations << "\ninitial_residual " << cg.initial_residual_norm
              << "\nresidual " << cg.residual_norm << "\nmin_K " << *cold << "\nmax_K " << *hot
              << "\ngenerated_W " << balance.generated_heat << "\nneumann_out_W " << balance.neumann_outward_heat
              << "\nconvection_out_W " << balance.robin_outward_heat
              << "\nfixed_out_W " << balance.constraint_outward_heat << "\nbalance_defect_W " << balance.signed_balance_defect
              << "\nmax_free_residual_W " << balance.max_free_equation_residual << '\n';
    return 0;
}
} // namespace

#ifdef _WIN32
int wmain(int argc, wchar_t** argv) {
#else
int main(int argc, char** argv) {
#endif
    try {
        const std::vector<std::filesystem::path> args(argv, argv + argc);
        const int status = run(args);
        std::cout.flush();
        if (!std::cout) throw std::runtime_error("Failed to write solver response.");
        return status;
    } catch (const std::exception& error) {
        std::cerr << "nova-run: " << error.what() << '\n';
        return 1;
    }
}
