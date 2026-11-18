#include "nova/fem/nonlinear_heat_2d.hpp"
#include "test_framework.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using nova::fem::H1TriangleOrder;
using nova::fem::NonlinearHeat2DProblem;
using nova::fem::PicardHeat2DOptions;

NonlinearHeat2DProblem make_nonlinear_problem(std::size_t subdivisions = 3U) {
    NonlinearHeat2DProblem problem;
    const auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, subdivisions, subdivisions
    );
    problem.thermal.mesh = generated.mesh;
    problem.thermal.volumetric_capacity =
        [](double, double) { return 1.0; };
    problem.thermal.source =
        [](double, double, double) { return 0.0; };
    problem.thermal.initial_temperature =
        [](double x, double y) { return 1.0 + x * x + y * y; };
    problem.conductivity = [](double, double, double temperature) {
        return 1.0 + temperature * temperature;
    };
    return problem;
}

} // namespace

NOVA_TEST(nonlinear_heat_p1_p2_preserve_constant_temperature) {
    auto problem = make_nonlinear_problem(2U);
    problem.thermal.initial_temperature =
        [](double, double) { return 4.25; };
    problem.conductivity = [](double, double, double temperature) {
        return 0.1 + temperature * temperature;
    };
    for (const H1TriangleOrder order : {
        H1TriangleOrder::P1, H1TriangleOrder::P2
    }) {
        PicardHeat2DOptions options;
        options.order = order;
        options.steps = 3U;
        options.time_step = 0.1;
        const auto result = nova::fem::solve_picard_heat_2d(problem, options);
        NOVA_REQUIRE(result.step_diagnostics.size() == options.steps);
        for (const double temperature : result.temperature) {
            NOVA_REQUIRE_NEAR(temperature, 4.25, 1.0e-8);
        }
        NOVA_REQUIRE_NEAR(result.final_time, 0.3, 1.0e-14);
    }
}

NOVA_TEST(nonlinear_heat_p1_p2_linear_spacetime_manufactured_solution) {
    constexpr double alpha = 0.2;
    for (const H1TriangleOrder order : {
        H1TriangleOrder::P1, H1TriangleOrder::P2
    }) {
        auto problem = make_nonlinear_problem(3U);
        problem.thermal.initial_temperature = [](double x, double y) {
            return x + y;
        };
        problem.thermal.source = [](double, double, double) {
            // T=x+y+t; div((1+alpha*T)*grad(T))=2*alpha.
            return 1.0 - 2.0 * alpha;
        };
        problem.conductivity = [](double, double, double temperature) {
            return 1.0 + alpha * temperature;
        };
        const auto generated = nova::fem::make_rectangle_mesh_2d(
            0.0, 1.0, 0.0, 1.0, 3U, 3U
        );
        for (const auto& edge : generated.boundary_edges) {
            problem.thermal.dirichlet_edges.push_back({
                edge,
                [](double x, double y, double time) {
                    return x + y + time;
                }
            });
        }

        PicardHeat2DOptions options;
        options.order = order;
        options.steps = 4U;
        options.time_step = 0.05;
        const auto result = nova::fem::solve_picard_heat_2d(problem, options);
        const nova::fem::H1FunctionSpace2D space(problem.thermal.mesh, order);
        const auto exact = space.interpolate([](double x, double y) {
            return x + y + 0.2;
        });
        NOVA_REQUIRE(result.temperature.size() == exact.size());
        for (std::size_t i = 0; i < exact.size(); ++i) {
            NOVA_REQUIRE_NEAR(result.temperature[i], exact[i], 2.0e-7);
        }
        NOVA_REQUIRE(result.linear_solvers.size() >= options.steps);
        for (const auto& step : result.step_diagnostics) {
            NOVA_REQUIRE(step.normalized_update <= 1.0);
        }
    }
}

NOVA_TEST(nonlinear_heat_picard_handles_nonuniform_temperature_field) {
    auto problem = make_nonlinear_problem(3U);
    PicardHeat2DOptions options;
    options.order = H1TriangleOrder::P2;
    options.steps = 1U;
    options.time_step = 0.02;
    options.max_iterations = 60U;
    options.absolute_tolerance = 1.0e-9;
    options.relative_tolerance = 1.0e-7;

    const nova::fem::H1FunctionSpace2D space(
        problem.thermal.mesh, options.order
    );
    const auto initial = space.interpolate(problem.thermal.initial_temperature);
    const auto result = nova::fem::solve_picard_heat_2d(problem, options);
    NOVA_REQUIRE(result.step_diagnostics.front().iterations > 1U);
    NOVA_REQUIRE(result.step_diagnostics.front().normalized_update <= 1.0);
    bool moved = false;
    for (std::size_t i = 0; i < initial.size(); ++i) {
        NOVA_REQUIRE(std::isfinite(result.temperature[i]));
        if (std::abs(result.temperature[i] - initial[i]) > 1.0e-4) {
            moved = true;
        }
    }
    NOVA_REQUIRE(moved);
}

NOVA_TEST(nonlinear_heat_rejects_nonpositive_conductivity) {
    auto problem = make_nonlinear_problem();
    problem.conductivity = [](double, double, double) { return -1.0; };
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_picard_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(nonlinear_heat_rejects_missing_material_or_bad_controls) {
    auto problem = make_nonlinear_problem();
    problem.conductivity = {};
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_picard_heat_2d(problem));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);

    problem = make_nonlinear_problem();
    PicardHeat2DOptions options;
    options.relaxation = 0.0;
    caught = false;
    try {
        static_cast<void>(nova::fem::solve_picard_heat_2d(problem, options));
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}

NOVA_TEST(nonlinear_heat_picard_iteration_budget_is_enforced) {
    auto problem = make_nonlinear_problem();
    PicardHeat2DOptions options;
    options.max_iterations = 1U;
    options.time_step = 0.2;
    options.absolute_tolerance = 1.0e-13;
    options.relative_tolerance = 0.0;
    bool caught = false;
    try {
        static_cast<void>(nova::fem::solve_picard_heat_2d(problem, options));
    } catch (const std::runtime_error&) {
        caught = true;
    }
    NOVA_REQUIRE(caught);
}
