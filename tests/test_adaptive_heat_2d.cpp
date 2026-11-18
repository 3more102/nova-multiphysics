#include "test_framework.hpp"

#include "nova/fem/adaptive_heat_2d.hpp"

#include <cmath>
#include <functional>
#include <map>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using nova::fem::AdaptiveHeat2DOptions;
using nova::fem::H1TriangleOrder;
using nova::fem::TransientHeat2DProblem;
using nova::fem::TransientTimeScheme;

TransientHeat2DProblem base_problem(std::size_t subdivisions = 2U) {
    auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, subdivisions, subdivisions
    );
    TransientHeat2DProblem problem;
    problem.mesh = std::move(generated.mesh);
    problem.volumetric_capacity = [](double, double) { return 1.0; };
    problem.conductivity = [](double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    problem.initial_temperature = [](double, double) { return 0.0; };
    return problem;
}

void set_boundary_trace(
    TransientHeat2DProblem& problem,
    const nova::fem::TransientScalarField2D& prescribed
) {
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> incidence;
    for (const auto& triangle : problem.mesh.elements) {
        ++incidence[std::minmax(triangle[0], triangle[1])];
        ++incidence[std::minmax(triangle[1], triangle[2])];
        ++incidence[std::minmax(triangle[2], triangle[0])];
    }
    for (const auto& [edge, number] : incidence) {
        if (number == 1U) {
            problem.dirichlet_edges.push_back({
                {edge.first, edge.second}, prescribed
            });
        }
    }
}

TransientHeat2DProblem cubic_manufactured_problem() {
    auto problem = base_problem(3U);
    // Exact P2 space/time solution u=t^3*x*(1-x).
    problem.source = [](double x, double, double t) {
        return 3.0 * t * t * x * (1.0 - x) + 2.0 * t * t * t;
    };
    set_boundary_trace(problem, [](double x, double, double t) {
        return t * t * t * x * (1.0 - x);
    });
    return problem;
}

double final_rms_error(
    const TransientHeat2DProblem& problem,
    const nova::fem::AdaptiveHeat2DResult& result,
    H1TriangleOrder order,
    const nova::fem::SpatialScalarField2D& expected
) {
    nova::fem::H1FunctionSpace2D space(problem.mesh, order);
    const auto target = space.interpolate(expected);
    NOVA_REQUIRE(target.size() == result.temperature.size());
    double sum = 0.0;
    for (std::size_t i = 0; i < target.size(); ++i) {
        const double delta = result.temperature[i] - target[i];
        sum += delta * delta;
    }
    return std::sqrt(sum / static_cast<double>(target.size()));
}

} // namespace

NOVA_TEST(transient_heat_start_time_and_exact_dof_restart_p1_p2) {
    for (const auto order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        auto problem = base_problem();
        nova::fem::H1FunctionSpace2D space(problem.mesh, order);
        problem.initial_dof_values = space.interpolate(
            [](double x, double y) { return x + y + 2.0; }
        );
        problem.initial_temperature = {};
        problem.source = [](double, double, double) { return 1.0; };
        set_boundary_trace(problem, [](double x, double y, double t) {
            return x + y + t;
        });
        nova::fem::TransientHeat2DOptions opts;
        opts.order = order;
        opts.scheme = TransientTimeScheme::CrankNicolson;
        opts.start_time = 2.0;
        opts.time_step = 0.2;
        opts.steps = 3U;
        const auto result = nova::fem::solve_transient_heat_2d(problem, opts);
        NOVA_REQUIRE_NEAR(result.final_time, 2.6, 1.0e-12);
        NOVA_REQUIRE(result.diagnostics.size() == 3U);
        NOVA_REQUIRE_NEAR(result.diagnostics.front().time, 2.2, 1.0e-12);
        const auto expected = space.interpolate(
            [](double x, double y) { return x + y + 2.6; }
        );
        for (std::size_t i = 0; i < expected.size(); ++i) {
            NOVA_REQUIRE_NEAR(result.temperature[i], expected[i], 2.0e-9);
        }
    }
}

NOVA_TEST(adaptive_heat_preserves_constant_p1_p2_with_zero_flux) {
    for (const auto order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        auto problem = base_problem();
        problem.initial_temperature = [](double, double) { return 3.25; };
        AdaptiveHeat2DOptions opts;
        opts.order = order;
        opts.end_time = 0.75;
        opts.initial_step = 0.07;
        opts.max_step = 0.25;
        opts.absolute_tolerance = 1.0e-7;
        const auto result = nova::fem::solve_adaptive_heat_2d(problem, opts);
        NOVA_REQUIRE_NEAR(result.final_time, opts.end_time, 1.0e-14);
        NOVA_REQUIRE(result.accepted_steps >= 1U);
        NOVA_REQUIRE(result.rejected_steps == 0U);
        NOVA_REQUIRE(result.diagnostics.size() == result.accepted_steps);
        NOVA_REQUIRE(result.substep_solvers.size() == 2U * result.accepted_steps);
        for (const double value : result.temperature) {
            NOVA_REQUIRE_NEAR(value, 3.25, 1.0e-9);
        }
        for (const auto& d : result.diagnostics) {
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 2.0e-9);
        }
    }
}

NOVA_TEST(adaptive_heat_crank_nicolson_rejects_large_steps_and_refines) {
    const auto problem = cubic_manufactured_problem();
    AdaptiveHeat2DOptions coarse;
    coarse.order = H1TriangleOrder::P2;
    coarse.scheme = TransientTimeScheme::CrankNicolson;
    coarse.end_time = 1.0;
    coarse.initial_step = 0.75;
    coarse.max_step = 0.75;
    coarse.min_step = 1.0e-5;
    coarse.absolute_tolerance = 1.0e-3;
    coarse.relative_tolerance = 0.0;

    AdaptiveHeat2DOptions fine = coarse;
    fine.absolute_tolerance = 1.0e-6;

    const auto loose = nova::fem::solve_adaptive_heat_2d(problem, coarse);
    const auto strict = nova::fem::solve_adaptive_heat_2d(problem, fine);

    NOVA_REQUIRE_NEAR(loose.final_time, 1.0, 1.0e-14);
    NOVA_REQUIRE_NEAR(strict.final_time, 1.0, 1.0e-14);
    NOVA_REQUIRE(strict.rejected_steps > 0U);
    NOVA_REQUIRE(strict.accepted_steps > loose.accepted_steps);
    NOVA_REQUIRE(strict.attempts.size() ==
        strict.accepted_steps + strict.rejected_steps);
    double previous_time = -1.0;
    for (const auto& attempt : strict.attempts) {
        NOVA_REQUIRE(attempt.start_time >= previous_time);
        previous_time = attempt.start_time;
        if (attempt.accepted) {
            NOVA_REQUIRE(attempt.normalized_error <= 1.0);
        } else {
            NOVA_REQUIRE(attempt.normalized_error > 1.0);
        }
    }
    for (const auto& d : strict.diagnostics) {
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 2.0e-8);
    }
    const auto exact = [](double x, double) { return x * (1.0 - x); };
    const double loose_error = final_rms_error(
        problem, loose, H1TriangleOrder::P2, exact
    );
    const double strict_error = final_rms_error(
        problem, strict, H1TriangleOrder::P2, exact
    );
    NOVA_REQUIRE(loose_error > strict_error);
    NOVA_REQUIRE(strict_error < 5.0e-5);
}

NOVA_TEST(adaptive_heat_backward_euler_reuses_operators_across_trials_p1_p2) {
    for (const auto order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        auto problem = base_problem(2U);
        std::size_t capacity_evaluations = 0U;
        std::size_t conductivity_evaluations = 0U;
        problem.volumetric_capacity = [&](double, double) {
            ++capacity_evaluations;
            return 1.0;
        };
        problem.conductivity = [&](double, double) {
            ++conductivity_evaluations;
            return 1.0;
        };
        // Exact uniform solution u=t^2 (unit volumetric capacity, zero grad).
        problem.source = [](double, double, double t) { return 2.0 * t; };
        AdaptiveHeat2DOptions opts;
        opts.order = order;
        opts.scheme = TransientTimeScheme::BackwardEuler;
        opts.end_time = 0.8;
        opts.initial_step = 0.4;
        opts.max_step = 0.4;
        opts.min_step = 1.0e-5;
        opts.absolute_tolerance = 2.0e-4;
        opts.relative_tolerance = 0.0;
        const auto result = nova::fem::solve_adaptive_heat_2d(problem, opts);
        NOVA_REQUIRE(result.accepted_steps > 2U);
        NOVA_REQUIRE(result.rejected_steps > 0U);
        NOVA_REQUIRE_NEAR(result.final_time, 0.8, 1.0e-14);
        for (const double temperature : result.temperature) {
            NOVA_REQUIRE_NEAR(temperature, 0.64, 0.03);
        }
        for (const auto& d : result.diagnostics) {
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 2.0e-9);
        }

        // Six triangle quadrature evaluations per cell and material field:
        // exact equality proves that M and K were assembled only once.
        const std::size_t expected_evaluations =
            6U * problem.mesh.elements.size();
        NOVA_REQUIRE(capacity_evaluations == expected_evaluations);
        NOVA_REQUIRE(conductivity_evaluations == expected_evaluations);
    }
}

NOVA_TEST(adaptive_heat_rejects_invalid_tolerance_and_time_window) {
    auto problem = base_problem(1U);
    AdaptiveHeat2DOptions opts;
    opts.end_time = opts.start_time;
    bool bad_interval = false;
    try {
        static_cast<void>(nova::fem::solve_adaptive_heat_2d(problem, opts));
    } catch (const std::invalid_argument&) {
        bad_interval = true;
    }
    NOVA_REQUIRE(bad_interval);
    opts.end_time = 1.0;
    opts.absolute_tolerance = 0.0;
    bool bad_tolerance = false;
    try {
        static_cast<void>(nova::fem::solve_adaptive_heat_2d(problem, opts));
    } catch (const std::invalid_argument&) {
        bad_tolerance = true;
    }
    NOVA_REQUIRE(bad_tolerance);
}

NOVA_TEST(adaptive_heat_enforces_maximum_trial_attempts) {
    auto problem = cubic_manufactured_problem();
    AdaptiveHeat2DOptions opts;
    opts.order = H1TriangleOrder::P2;
    opts.end_time = 1.0;
    opts.initial_step = 0.75;
    opts.max_step = 0.75;
    opts.absolute_tolerance = 1.0e-10;
    opts.relative_tolerance = 0.0;
    opts.max_attempts = 1U;
    bool exhausted = false;
    try {
        static_cast<void>(nova::fem::solve_adaptive_heat_2d(problem, opts));
    } catch (const std::runtime_error&) {
        exhausted = true;
    }
    NOVA_REQUIRE(exhausted);
}

NOVA_TEST(adaptive_heat_rejects_nonfinite_restart_state) {
    auto problem = base_problem(1U);
    nova::fem::H1FunctionSpace2D space(problem.mesh, H1TriangleOrder::P1);
    problem.initial_dof_values = space.interpolate([](double, double) { return 1.0; });
    problem.initial_dof_values[0] = std::numeric_limits<double>::infinity();
    bool rejected = false;
    try {
        static_cast<void>(nova::fem::solve_adaptive_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}
