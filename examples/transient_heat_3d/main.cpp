#include "nova/fem/transient_heat_3d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
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
