#include "nova/fem/transient_heat_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/io/vtk_series.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

bool safe_series_prefix(std::string_view name) {
    if (name.empty() || name.size() > 80U) return false;
    for (const char c : name) {
        const bool letter = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
        const bool number = c >= '0' && c <= '9';
        if (!(letter || number || c == '_' || c == '-')) return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 1 && (argc != 3 || std::string_view(argv[1]) != "--vtk-series" ||
                      !safe_series_prefix(argv[2]))) {
        std::cerr << "Usage: nova_transient_heat_3d [--vtk-series safe_prefix]\\n";
        return 2;
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
    const auto result = solve_transient_heat_3d(problem, options);
    if (argc == 3) {
        const std::string prefix(argv[2]);
        std::vector<nova::io::VtkTimeSeriesEntry> frames;
        frames.reserve(options.steps + 1U);
        const auto write_frame = [&](std::size_t index, double time,
                                     const std::vector<double>& values) {
            std::ostringstream name;
            name << prefix << '_' << std::setfill('0') << std::setw(6) << index << ".vtk";
            std::ofstream output(name.str(), std::ios::out | std::ios::trunc);
            if (!output) throw std::runtime_error("Cannot create 3D VTK frame.");
            nova::io::write_legacy_vtk_scalar_3d(output, problem.diffusion.mesh, values);
            output.close();
            if (!output) throw std::runtime_error("Failed to save 3D VTK frame.");
            frames.push_back({time, name.str()});
        };
        // Example initial field is already consistent with its fixed 300 K boundary.
        write_frame(0U, options.start_time, problem.initial_temperature);
        for (std::size_t step = 0; step < result.linear_solvers.size(); ++step) {
            write_frame(step + 1U, result.diagnostics[step].time,
                        result.linear_solvers[step].solution);
        }
        const std::string collection = prefix + ".pvd";
        std::ofstream manifest(collection, std::ios::out | std::ios::trunc);
        if (!manifest) throw std::runtime_error("Cannot create PVD collection.");
        nova::io::write_pvd_collection(manifest, frames);
        manifest.close();
        if (!manifest) throw std::runtime_error("Failed to save PVD collection.");
        std::cout << "ParaView time series: " << collection << '\\n';
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
    return 0;
}
