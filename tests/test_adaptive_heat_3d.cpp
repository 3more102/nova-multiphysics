#include "test_framework.hpp"
#include "nova/fem/adaptive_heat_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
using namespace nova::fem;

TransientHeat3DProblem tetra_problem() {
    TransientHeat3DProblem p;
    p.diffusion.mesh = {{{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{{0,1,2,3}}};
    p.diffusion.cell_diffusivity = {1.0};
    p.cell_volumetric_capacity = {1.0};
    p.initial_temperature = {0.0,0.0,0.0,0.0};
    return p;
}

template<class F> bool throws_invalid(F&& f) {
    try { f(); } catch (const std::invalid_argument&) { return true; }
    return false;
}
} // namespace

NOVA_TEST(adaptive_heat_3d_insulated_uniform_source_and_energy) {
    auto p = tetra_problem();
    p.cell_volumetric_capacity = {4.0};
    p.diffusion.source = [](double, double, double) { return 8.0; };
    p.initial_temperature = {300.0, 300.0, 300.0, 300.0};
    AdaptiveHeat3DOptions options;
    options.scheme = TransientHeat3DScheme::CrankNicolson;
    options.end_time = 0.5;
    options.initial_step = 0.2;
    options.max_step = 0.2;
    options.absolute_tolerance = 1.0e-5;
    options.relative_tolerance = 0.0;
    const auto result = solve_adaptive_heat_3d(p, options);
    NOVA_REQUIRE_NEAR(result.final_time, 0.5, 1e-14);
    NOVA_REQUIRE(result.accepted_steps > 0);
    NOVA_REQUIRE(result.rejected_steps == 0);
    NOVA_REQUIRE(result.diagnostics.size() == result.accepted_steps);
    NOVA_REQUIRE(result.substep_solvers.size() == 2 * result.accepted_steps);
    for (double t : result.temperature) NOVA_REQUIRE_NEAR(t, 301.0, 1e-8);
    for (const auto& diag : result.diagnostics) {
        NOVA_REQUIRE_NEAR(diag.balance_residual, 0.0, 1e-9);
    }
}

NOVA_TEST(adaptive_heat_3d_rejections_and_decay_against_analytical_mode) {
    auto p = tetra_problem();
    p.diffusion.dirichlet_constraints = {{0,0},{1,0},{2,0}};
    p.initial_temperature[3] = 10.0;
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                               TransientHeat3DScheme::CrankNicolson}) {
        AdaptiveHeat3DOptions options;
        options.scheme = scheme;
        options.start_time = 3.0;
        options.end_time = 3.2;
        options.initial_step = 0.2;
        options.max_step = 0.2;
        options.min_step = 1e-5;
        options.absolute_tolerance = 2e-3;
        options.relative_tolerance = 0.0;
        options.solver_options.relative_tolerance = 1e-13;
        const auto r = solve_adaptive_heat_3d(p, options);
        NOVA_REQUIRE(r.rejected_steps > 0);
        NOVA_REQUIRE(r.accepted_steps > 0);
        NOVA_REQUIRE_NEAR(r.final_time, 3.2, 2e-14);
        NOVA_REQUIRE(r.attempts.size() == r.accepted_steps + r.rejected_steps);
        NOVA_REQUIRE(r.diagnostics.size() == r.accepted_steps);
        for (const auto& attempt : r.attempts) {
            NOVA_REQUIRE(attempt.accepted == (attempt.normalized_error <= 1.0));
        }
        for (const auto& d : r.diagnostics) NOVA_REQUIRE_NEAR(d.balance_residual, 0, 1e-9);
        for (std::size_t j=0; j<3; ++j) NOVA_REQUIRE_NEAR(r.temperature[j], 0.0, 1e-12);
        // Exact single unconstrained mode: T(3) = 10 exp(-10*0.2).
        NOVA_REQUIRE_NEAR(r.temperature[3], 10.0 * std::exp(-2.0), 0.08);
    }
}

NOVA_TEST(adaptive_heat_3d_restart_matches_fixed_step_when_no_rejections) {
    auto p = tetra_problem();
    p.diffusion.source = [](double, double, double) { return 2.0; };
    p.initial_temperature = {10,10,10,10};
    AdaptiveHeat3DOptions adaptive;
    adaptive.end_time = 0.3;
    adaptive.initial_step = 0.1;
    adaptive.max_step = 0.1;
    adaptive.absolute_tolerance = 1e-6;
    const auto a = solve_adaptive_heat_3d(p, adaptive);
    TransientHeat3DOptions fixed;
    fixed.scheme = adaptive.scheme;
    fixed.time_step = 0.05;
    fixed.steps = 6;
    const auto b = solve_transient_heat_3d(p, fixed);
    NOVA_REQUIRE(a.rejected_steps == 0);
    NOVA_REQUIRE(a.accepted_steps == 3);
    for (std::size_t i=0; i<4; ++i) NOVA_REQUIRE_NEAR(a.temperature[i], b.temperature[i], 1e-10);
}

NOVA_TEST(adaptive_heat_3d_rejects_invalid_options_and_minimum_step) {
    auto p = tetra_problem();
    AdaptiveHeat3DOptions opt;
    opt.end_time = 0.2;
    NOVA_REQUIRE(throws_invalid([&] {
        auto x=opt; x.end_time=0; (void)solve_adaptive_heat_3d(p,x);
    }));
    NOVA_REQUIRE(throws_invalid([&] {
        auto x=opt; x.min_step=0; (void)solve_adaptive_heat_3d(p,x);
    }));
    NOVA_REQUIRE(throws_invalid([&] {
        auto x=opt; x.absolute_tolerance=std::numeric_limits<double>::infinity();
        (void)solve_adaptive_heat_3d(p,x);
    }));
    NOVA_REQUIRE(throws_invalid([&] {
        auto x=opt; x.max_attempts=0; (void)solve_adaptive_heat_3d(p,x);
    }));
    NOVA_REQUIRE(throws_invalid([&] {
        auto x=opt; x.scheme=static_cast<TransientHeat3DScheme>(99);
        (void)solve_adaptive_heat_3d(p,x);
    }));
    p.diffusion.dirichlet_constraints = {{0,0},{1,0},{2,0}};
    p.initial_temperature[3] = 10.0;
    opt.scheme = TransientHeat3DScheme::BackwardEuler;
    opt.initial_step = 0.2;
    opt.min_step = 0.2;
    opt.max_step = 0.2;
    opt.absolute_tolerance = 1e-12;
    bool failed_min = false;
    try { (void)solve_adaptive_heat_3d(p,opt); }
    catch (const std::runtime_error&) { failed_min=true; }
    NOVA_REQUIRE(failed_min);
}
