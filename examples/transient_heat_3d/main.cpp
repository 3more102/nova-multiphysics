#include "nova/fem/transient_heat_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/io/vtk_series.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
    using namespace nova::fem;
    namespace fs = std::filesystem;

    const bool series = argc == 3 && std::string(argv[1]) == "--vtk-series";
    if (!series && argc != 1) {
        throw std::invalid_argument(
            "Usage: nova_transient_heat_3d [--vtk-series destination.pvd]");
    }
    const fs::path destination = series ? fs::path(argv[2]) : fs::path{};
    std::string prefix;
    if (series) {
        prefix = destination.stem().string();
        if (destination.extension() != ".pvd" || prefix.empty() ||
            prefix.front() == '.' || prefix.find("..") != std::string::npos ||
            !std::all_of(prefix.begin(), prefix.end(), [](char c) {
                return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                       (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.';
            })) {
            throw std::invalid_argument("Use a safe .pvd destination basename.");
        }
    }

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
    std::vector<nova::io::VtkTimeSeriesEntry> frames;
    if (series) {
        frames.reserve(options.steps + 1U);
        options.retain_step_solutions = false;
        options.on_step = [&](std::size_t step, double time,
                              const std::vector<double>& temperature) {
            std::ostringstream frame_name;
            frame_name << prefix << '_' << std::setw(6) << std::setfill('0')
                       << step << ".vtk";
            const auto name = frame_name.str();
            const fs::path file = destination.parent_path() / name;
            std::ofstream output(file);
            if (!output) throw std::runtime_error("Cannot open VTK frame: " + file.string());
            nova::io::write_legacy_vtk_scalar_3d(
                output, problem.diffusion.mesh, temperature, "temperature");
            output.close();
            if (!output) throw std::runtime_error("Cannot finish VTK frame: " + file.string());
            frames.push_back({time, name});
        };
    }
    const auto result = solve_transient_heat_3d(problem, options);
    if (series) {
        // Publish the collection only if all time-step solves and frames succeeded.
        std::ofstream output(destination);
        if (!output) throw std::runtime_error("Cannot open PVD collection.");
        nova::io::write_pvd_collection(output, frames);
        output.close();
        if (!output) throw std::runtime_error("Cannot finish PVD collection.");
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
              << "Final-step CG iterations: " << result.linear_solvers.back().iterations << '\\n';
    if (series) {
        std::cout << "ParaView collection: " << destination.string()
                  << " (" << frames.size() << " frames)\\n";
    }
    return 0;
    } catch (const std::exception& error) {
        std::cerr << "NOVA transient heat error: " << error.what() << '\\n';
        return 1;
    }
}
