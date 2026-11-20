#include "nova/fem/adaptive_heat_3d.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
    using namespace nova::fem;

    TransientHeat3DProblem problem;
    problem.diffusion.mesh = {
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
        {{0, 1, 2, 3}}
    };
    problem.diffusion.cell_diffusivity = {1.0};
    problem.diffusion.source = [](double, double, double) { return 5.0; };
    problem.cell_volumetric_capacity = {1.0};
    problem.initial_temperature = {300.0, 300.0, 300.0, 300.0};
    problem.diffusion.dirichlet_constraints = {{0, 300.0}};

    AdaptiveHeat3DOptions options;
    options.scheme = TransientHeat3DScheme::CrankNicolson;
    options.end_time = 0.5;
    options.initial_step = 0.25;
    options.min_step = 1.0e-5;
    options.max_step = 0.25;
    options.absolute_tolerance = 1.0e-5;
    options.relative_tolerance = 1.0e-4;
    options.solver_options.preconditioner =
        nova::linalg::PreconditionerType::IncompleteCholesky0;
    options.solver_options.relative_tolerance = 1.0e-12;

    std::size_t emitted_fields = 0U;
    options.on_step = [&](std::size_t, double, const std::vector<double>&) {
        ++emitted_fields;
    };

    const auto result = solve_adaptive_heat_3d(problem, options);
    const auto [minimum, maximum] = std::minmax_element(
        result.temperature.begin(), result.temperature.end());
    double largest_error = 0.0;
    for (const auto& attempt : result.attempts) {
        largest_error = std::max(largest_error, attempt.normalized_error);
    }
    double largest_balance = 0.0;
    for (const auto& diagnostic : result.diagnostics) {
        largest_balance = std::max(largest_balance,
                                   std::abs(diagnostic.balance_residual));
    }
    for (const auto& solver : result.substep_solvers) {
        if (!solver.converged || !solver.solution.empty()) return 1;
    }

    std::cout << std::setprecision(10)
              << "NOVA adaptive 3D tetrahedral heat demonstration\n"
              << "final_time " << result.final_time << '\n'
              << "accepted_steps " << result.accepted_steps << '\n'
              << "rejected_steps " << result.rejected_steps << '\n'
              << "emitted_fields " << emitted_fields << '\n'
              << "temperature_range " << *minimum << " " << *maximum << '\n'
              << "largest_normalized_attempt_error " << largest_error << '\n'
              << "max_abs_energy_balance_residual " << largest_balance << '\n'
              << "full_step_ic0_factor_builds "
              << result.full_step_ic0_factor_builds << '\n'
              << "half_step_ic0_factor_builds "
              << result.half_step_ic0_factor_builds << '\n';
    return std::abs(result.final_time - options.end_time) <= 1.0e-12 ? 0 : 2;
}
