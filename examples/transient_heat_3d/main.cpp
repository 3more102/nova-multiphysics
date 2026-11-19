#include "nova/fem/transient_heat_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/io/vtk_series.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

int main(int argc, char** argv) {
    using namespace nova::fem;
    const bool vtk_series = argc == 3 && std::string_view(argv[1]) == "--vtk-series";
    if (argc != 1 && !vtk_series) {
        std::cerr << "Usage: nova_transient_heat_3d [--vtk-series safe_prefix]\n";
        return 2;
    }
    std::string prefix;
    if (vtk_series) {
        prefix = argv[2];
        if (prefix.empty() || prefix.front() == '.' ||
            prefix.size() > 100U ||
            !std::all_of(prefix.begin(), prefix.end(), [](char c) {
                return (c >= 'a' && c <= 'z') ||
                       (c >= 'A' && c <= 'Z') ||
                       (c >= '0' && c <= '9') || c == '_' || c == '-';
            })) {
            std::cerr << "PVD prefix must be a safe ASCII basename (maximum 100 characters).\n";
            return 2;
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
    if (vtk_series) {
        options.state_observer = [&](std::size_t step, double time,
                                      const std::vector<double>& temperatures) {
            std::ostringstream name_builder;
            name_builder << prefix << '_' << std::setfill('0') << std::setw(6)
                         << step << ".vtk";
            const std::string frame_name = name_builder.str();
            std::ofstream stream(frame_name, std::ios::binary);
            if (!stream) throw std::runtime_error("Cannot open thermal VTK frame.");
            nova::io::write_legacy_vtk_scalar_3d(
                stream, problem.diffusion.mesh, temperatures);
            stream.close();
            if (!stream) throw std::runtime_error("Cannot finish thermal VTK frame.");
            frames.push_back({time, frame_name});
        };
    }
    const auto result = solve_transient_heat_3d(problem, options);
    if (vtk_series) {
        std::ofstream collection(prefix + ".pvd", std::ios::binary);
        if (!collection) throw std::runtime_error("Cannot open thermal PVD collection.");
        nova::io::write_pvd_collection(collection, frames);
        collection.close();
        if (!collection) throw std::runtime_error("Cannot finish thermal PVD collection.");
    }
    const auto [lo, hi] = std::minmax_element(
        result.temperature.begin(), result.temperature.end()
    );
    const auto& energy = result.diagnostics.back();
    if (vtk_series) {
        std::cout << "PVD collection written: " << prefix << ".pvd ("
                  << frames.size() << " frames)\n";
    }
    std::cout << std::setprecision(12)
              << "NOVA 3D transient heat demonstration\n"
              << "P1 tetrahedra: " << problem.diffusion.mesh.elements.size() << '\n'
              << "Time [s]: " << result.final_time << '\n'
              << "Temperature range [K]: " << *lo << " ... " << *hi << '\n'
              << "Final-step energy residual [J]: " << energy.balance_residual << '\n'
              << "Final-step CG iterations: " << result.linear_solvers.back().iterations << '\n';
    return 0;
}
