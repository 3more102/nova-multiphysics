#include "nova/fem/nonlinear_heat_2d.hpp"

#include <algorithm>
#include <exception>
#include <iostream>

int main() {
    try {
        nova::fem::NonlinearHeat2DProblem problem;
        problem.thermal.mesh = nova::fem::make_rectangle_mesh_2d(
            0.0, 1.0, 0.0, 1.0, 5U, 5U
        ).mesh;
        problem.thermal.volumetric_capacity =
            [](double, double) { return 1.0; };
        problem.thermal.initial_temperature =
            [](double x, double y) {
                return 1.0 + 2.0 * x * (1.0 - x) * y * (1.0 - y);
            };
        problem.thermal.source =
            [](double, double, double) { return 0.1; };
        // All edges are insulated by default.
        problem.conductivity =
            [](double, double, double temperature) {
                return 1.0 + 0.3 * temperature * temperature;
            };

        nova::fem::PicardHeat2DOptions options;
        options.order = nova::fem::H1TriangleOrder::P2;
        options.time_step = 0.05;
        options.steps = 4U;
        const auto result = nova::fem::solve_picard_heat_2d(problem, options);
        const auto [lo, hi] = std::minmax_element(
            result.temperature.begin(), result.temperature.end()
        );
        std::cout << "NOVA nonlinear thermal Picard demo\n";
        std::cout << "t = " << result.final_time
                  << ", T_min = " << *lo << ", T_max = " << *hi << '\n';
        for (const auto& step : result.step_diagnostics) {
            std::cout << "t = " << step.time
                      << ", Picard iterations = " << step.iterations
                      << ", scaled update = " << step.normalized_update << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Nonlinear thermal solve failed: " << error.what() << '\n';
        return 1;
    }
}
