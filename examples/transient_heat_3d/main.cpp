#include "nova/fem/transient_heat_3d.hpp"
#include "nova/io/pvd_collection.hpp"
#include "nova/io/vtk_legacy_3d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

int main(int argc, char** argv) {
    try {
    const bool series = argc == 3 && std::string_view(argv[1]) == "--vtk-series";
    if (argc != 1 && !series) {
        throw std::invalid_argument("Usage: nova_transient_heat_3d [--vtk-series output.pvd]");
    }
    using namespace nova::fem;

    // A bonded 1-metre cube initially at 300 K, internally heated and
    // cooled at x=1. Other faces are insulated; x=0 is held at 300 K.
    auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 3, 3, 3);
    TransientHeat3DProblem problem;
    problem.diffusion.mesh = std::move(box.mesh);
    problem.diffusion.cell_diffusivity.assign(
        problem.diffusion.mesh.elements.size(), 15.0
    );
    problem.diffusion.source = [](double, double, double) { return 2.0e5; };
    problem.volumetric_capacity =
        [](double, double, double) { return 3.0e6; };
    problem.initial_temperature.assign(
        problem.diffusion.mesh.nodes.size(), 300.0
    );

    for (const auto& face : box.boundary_faces) {
        bool cold_wall = true;
        bool fixed_wall = true;
        for (std::size_t node : face) {
            const double x = problem.diffusion.mesh.nodes[node].x;
            cold_wall = cold_wall && x == 1.0;
            fixed_wall = fixed_wall && x == 0.0;
        }
        if (fixed_wall) {
            problem.diffusion.dirichlet_faces.push_back({
                face, [](double, double, double) { return 300.0; }
            });
        } else if (cold_wall) {
            problem.diffusion.robin_faces.push_back({
                face,
                [](double, double, double) { return 30.0; },
                [](double, double, double) { return 295.0; }
            });
        }
    }

    TransientHeat3DOptions options;
    options.time_step = 20.0;
    options.steps = 5;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::Jacobi;
    options.solver_options.relative_tolerance = 1e-11;

    std::filesystem::path collection_path;
    std::vector<nova::io::ParaviewTimeFrame> frames;
    if (series) {
        collection_path = argv[2];
        const std::string stem = collection_path.stem().string();
        const auto safe_stem = [](const std::string& name) {
            if (name.empty() || name == "." || name == "..") return false;
            return std::all_of(name.begin(), name.end(), [](unsigned char c) {
                return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                       (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.';
            });
        };
        if (collection_path.extension() != ".pvd" || !safe_stem(stem)) {
            throw std::invalid_argument("Expected an ASCII-safe .pvd output filename.");
        }
        options.retain_step_solutions = false;
        options.step_observer = [&, stem](double time, const std::vector<double>& values) {
            std::ostringstream suffix;
            suffix << std::setfill('0') << std::setw(6) << frames.size();
            const std::string name = stem + "_" + suffix.str() + ".vtk";
            std::ofstream frame(collection_path.parent_path() / name);
            if (!frame) throw std::runtime_error("Cannot create VTK time-series frame.");
            nova::io::write_legacy_vtk_scalar_3d(
                frame, problem.diffusion.mesh, values, "temperature");
            frame.close();
            if (!frame) throw std::runtime_error("Cannot close VTK time-series frame.");
            frames.push_back({time, name});
        };
    }
    const auto result = solve_transient_heat_3d(problem, options);
    if (series) {
        std::ofstream collection(collection_path);
        if (!collection) throw std::runtime_error("Cannot create ParaView collection.");
        nova::io::write_paraview_collection(collection, frames);
        collection.close();
        if (!collection) throw std::runtime_error("Cannot close ParaView collection.");
    }
    const auto [lo, hi] = std::minmax_element(
        result.temperature.begin(), result.temperature.end()
    );
    const auto& energy = result.diagnostics.back();
    std::cout << std::setprecision(12)
              << "NOVA 3D transient heat demonstration\n"
              << "P1 tetrahedra: " << problem.diffusion.mesh.elements.size() << '\n'
              << "Time [s]: " << result.final_time << '\n'
              << "Temperature range [K]: " << *lo << " ... " << *hi << '\n'
              << "Final-step energy residual [J]: " << energy.balance_residual << '\n'
              << "Final-step CG iterations: " << result.linear_solvers.back().iterations << '\n';
    if (series) std::cout << "ParaView frames: " << frames.size() << '\n';
    return 0;
    } catch (const std::exception& error) {
        std::cerr << "NOVA transient heat error: " << error.what() << '\n';
        return 1;
    }
}
