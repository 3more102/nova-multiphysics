#include "test_framework.hpp"
#include "nova/fem/transient_heat_3d.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
nova::fem::TransientHeat3DProblem make_streamed_unit_heat() {
    using namespace nova::fem;
    TransientHeat3DProblem p;
    p.diffusion.mesh = {{{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{{0,1,2,3}}};
    p.diffusion.cell_diffusivity = {1.0};
    p.diffusion.source = [](double,double,double) { return 8.0; };
    p.cell_volumetric_capacity = {4.0};
    p.initial_temperature = {300,300,300,300};
    return p;
}
} // namespace

NOVA_TEST(transient_heat_3d_streaming_observer_sees_exact_initial_and_every_step) {
    using namespace nova::fem;
    for (auto scheme : {TransientHeat3DScheme::BackwardEuler,
                        TransientHeat3DScheme::CrankNicolson}) {
        auto p = make_streamed_unit_heat();
        TransientHeat3DOptions options;
        options.scheme = scheme;
        options.time_step = 0.125;
        options.start_time = 2.0;
        options.steps = 4;
        options.retain_step_solutions = false;
        options.solver_options.relative_tolerance = 1e-13;
        options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
        std::size_t seen = 0;
        options.on_step = [&](std::size_t step, double time,
                              const std::vector<double>& field) {
            NOVA_REQUIRE(step == seen++);
            NOVA_REQUIRE_NEAR(time, 2.0 + 0.125 * static_cast<double>(step), 1e-14);
            NOVA_REQUIRE(field.size() == 4U);
            for (double t : field) {
                NOVA_REQUIRE_NEAR(t, 300.0 + 0.25 * static_cast<double>(step), 1e-9);
            }
        };
        const auto out = solve_transient_heat_3d(p, options);
        NOVA_REQUIRE(seen == options.steps + 1U);
        NOVA_REQUIRE_NEAR(out.final_time, 2.5, 1e-14);
        NOVA_REQUIRE(out.linear_solvers.size() == options.steps);
        NOVA_REQUIRE(out.diagnostics.size() == options.steps);
        for (const auto& solver : out.linear_solvers) {
            NOVA_REQUIRE(solver.converged);
            NOVA_REQUIRE(solver.solution.empty());
        }
        for (double temperature : out.temperature) {
            NOVA_REQUIRE_NEAR(temperature, 301.0, 1e-9);
        }
        for (const auto& d : out.diagnostics) {
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1e-10);
        }
    }
}

NOVA_TEST(transient_heat_3d_streaming_respects_initial_dirichlet_projection_and_default_history) {
    using namespace nova::fem;
    auto p = make_streamed_unit_heat();
    p.diffusion.source = [](double,double,double) { return 0.0; };
    p.initial_temperature = {100,10,10,10};
    p.diffusion.dirichlet_constraints = {{0, 0.0}};
    TransientHeat3DOptions options;
    options.steps = 2;
    options.time_step = 0.1;
    options.solver_options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    std::vector<std::vector<double>> snapshots;
    options.on_step = [&](std::size_t, double, const std::vector<double>& field) {
        snapshots.push_back(field);
    };
    const auto out = solve_transient_heat_3d(p, options);
    NOVA_REQUIRE(snapshots.size() == 3U);
    NOVA_REQUIRE_NEAR(snapshots[0][0], 0.0, 1e-15);
    NOVA_REQUIRE_NEAR(snapshots[0][1], 10.0, 1e-15);
    NOVA_REQUIRE(out.linear_solvers.size() == 2U);
    NOVA_REQUIRE(out.linear_solvers[0].solution.size() == 4U);
    NOVA_REQUIRE(out.linear_solvers[1].solution.size() == 4U);
    for (std::size_t i = 0; i < 4U; ++i) {
        NOVA_REQUIRE_NEAR(snapshots[1][i], out.linear_solvers[0].solution[i], 1e-13);
        NOVA_REQUIRE_NEAR(snapshots[2][i], out.temperature[i], 1e-13);
        NOVA_REQUIRE_NEAR(snapshots[2][i], out.linear_solvers[1].solution[i], 1e-13);
    }
}

NOVA_TEST(transient_heat_3d_streaming_callback_exceptions_propagate) {
    using namespace nova::fem;
    auto p = make_streamed_unit_heat();
    TransientHeat3DOptions options;
    options.steps = 2;
    bool initial_threw = false;
    options.on_step = [](std::size_t, double, const std::vector<double>&) {
        throw std::runtime_error("initial observer exception");
    };
    try {
        (void)solve_transient_heat_3d(p, options);
    } catch (const std::runtime_error& err) {
        initial_threw = std::string(err.what()) == "initial observer exception";
    }
    NOVA_REQUIRE(initial_threw);

    bool step_threw = false;
    options.on_step = [](std::size_t step, double, const std::vector<double>&) {
        if (step == 1U) throw std::runtime_error("step observer exception");
    };
    try {
        (void)solve_transient_heat_3d(p, options);
    } catch (const std::runtime_error& err) {
        step_threw = std::string(err.what()) == "step observer exception";
    }
    NOVA_REQUIRE(step_threw);

    std::size_t calls = 0;
    options.time_step = -1.0;
    options.on_step = [&](std::size_t, double, const std::vector<double>&) {
        ++calls;
    };
    try {
        (void)solve_transient_heat_3d(p, options);
    } catch (const std::invalid_argument&) {}
    NOVA_REQUIRE(calls == 0U);
}


NOVA_TEST(transient_heat_3d_rejects_time_steps_that_do_not_advance_physical_time) {
    using namespace nova::fem;
    auto p = make_streamed_unit_heat();
    TransientHeat3DOptions options;
    options.start_time = 1.0;
    options.time_step = 0.5 * std::numeric_limits<double>::epsilon();
    options.steps = 2;
    std::size_t callbacks = 0;
    options.on_step = [&](std::size_t, double, const std::vector<double>&) {
        ++callbacks;
    };

    bool rejected = false;
    try {
        (void)solve_transient_heat_3d(p, options);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
    NOVA_REQUIRE(callbacks == 0U);
}
