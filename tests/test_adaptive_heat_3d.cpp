#include "test_framework.hpp"
#include "nova/fem/adaptive_heat_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace nova::fem;

TransientHeat3DProblem base_tetra() {
    TransientHeat3DProblem p;
    p.diffusion.mesh = {
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
        {{0, 1, 2, 3}}
    };
    p.diffusion.cell_diffusivity = {1.0};
    p.cell_volumetric_capacity = {1.0};
    p.initial_temperature = {0.0, 0.0, 0.0, 0.0};
    return p;
}
} // namespace

NOVA_TEST(adaptive_heat_3d_exact_uniform_heating_and_accepted_observer) {
    for (const auto scheme : {TransientHeat3DScheme::BackwardEuler,
                              TransientHeat3DScheme::CrankNicolson}) {
        auto p = base_tetra();
        p.cell_volumetric_capacity = {4.0};
        p.diffusion.cell_source = {8.0};
        p.initial_temperature = {300.0, 300.0, 300.0, 300.0};
        AdaptiveHeat3DOptions o;
        o.scheme = scheme;
        o.start_time = 0.3;
        o.end_time = 1.1;
        o.initial_step = 0.17;
        o.max_step = 0.35;
        std::vector<double> observed_times;
        std::vector<std::vector<double>> observed_states;
        o.on_step = [&](std::size_t index, double t, const std::vector<double>& values) {
            NOVA_REQUIRE(index == observed_times.size());
            observed_times.push_back(t);
            observed_states.push_back(values);
        };
        const auto r = solve_adaptive_heat_3d(p, o);
        NOVA_REQUIRE(r.accepted_steps >= 3U);
        NOVA_REQUIRE(r.rejected_steps == 0U);
        NOVA_REQUIRE(r.attempts.size() == r.accepted_steps);
        NOVA_REQUIRE(r.diagnostics.size() == r.accepted_steps);
        NOVA_REQUIRE(observed_times.size() == r.accepted_steps + 1U);
        NOVA_REQUIRE_NEAR(observed_times.front(), 0.3, 1e-12);
        NOVA_REQUIRE_NEAR(r.final_time, 1.1, 1e-12);
        NOVA_REQUIRE_NEAR(observed_times.back(), r.final_time, 1e-12);
        for (std::size_t i = 0; i < observed_times.size(); ++i) {
            if (i != 0U) NOVA_REQUIRE(observed_times[i] > observed_times[i - 1U]);
            const double expected = 300.0 + 2.0 * (observed_times[i] - 0.3);
            for (double v : observed_states[i]) NOVA_REQUIRE_NEAR(v, expected, 2e-8);
        }
        for (double v : r.temperature) NOVA_REQUIRE_NEAR(v, 301.6, 2e-8);
        for (const auto& d : r.diagnostics) {
            NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 2e-9);
        }

        // Exact fixed-mesh nodal restart from an accepted field.
        auto restarted = p;
        restarted.initial_temperature = r.temperature;
        o.on_step = {};
        o.start_time = r.final_time;
        o.end_time = 1.5;
        const auto continued = solve_adaptive_heat_3d(restarted, o);
        for (double v : continued.temperature) NOVA_REQUIRE_NEAR(v, 302.4, 2e-8);
    }
}

NOVA_TEST(adaptive_heat_3d_rejects_stiff_decay_and_tightening_improves_accuracy) {
    auto p = base_tetra();
    p.initial_temperature = {0.0, 0.0, 0.0, 10.0};
    p.diffusion.dirichlet_constraints = {{0, 0.0}, {1, 0.0}, {2, 0.0}};
    AdaptiveHeat3DOptions loose;
    loose.scheme = TransientHeat3DScheme::BackwardEuler;
    loose.end_time = 0.5;
    loose.initial_step = 0.5;
    loose.max_step = 0.5;
    loose.min_step = 1e-6;
    loose.absolute_tolerance = 0.05;
    loose.relative_tolerance = 0.0;
    std::size_t observed = 0;
    loose.on_step = [&](std::size_t index, double, const std::vector<double>&) {
        NOVA_REQUIRE(index == observed++);
    };
    const auto a = solve_adaptive_heat_3d(p, loose);
    NOVA_REQUIRE(a.rejected_steps > 0U);
    NOVA_REQUIRE(observed == a.accepted_steps + 1U);
    auto strict = loose;
    strict.absolute_tolerance = 0.001;
    strict.on_step = {};
    const auto b = solve_adaptive_heat_3d(p, strict);
    NOVA_REQUIRE(b.accepted_steps > a.accepted_steps);
    NOVA_REQUIRE(b.rejected_steps > 0U);
    const double exact = 10.0 * std::exp(-5.0);
    NOVA_REQUIRE(std::abs(b.temperature[3] - exact) <
                 std::abs(a.temperature[3] - exact));
    for (const auto& d : b.diagnostics) {
        NOVA_REQUIRE_NEAR(d.balance_residual, 0.0, 1e-9);
    }
}

NOVA_TEST(adaptive_heat_3d_rejects_bad_controls_and_minimum_step_failure) {
    const auto p = base_tetra();
    AdaptiveHeat3DOptions o;
    o.end_time = o.start_time;
    bool bad_time = false;
    try { (void)solve_adaptive_heat_3d(p, o); }
    catch (const std::invalid_argument&) { bad_time = true; }
    NOVA_REQUIRE(bad_time);

    o.end_time = 0.5;
    o.absolute_tolerance = 0.0;
    o.relative_tolerance = 0.0;
    bool bad_tolerance = false;
    try { (void)solve_adaptive_heat_3d(p, o); }
    catch (const std::invalid_argument&) { bad_tolerance = true; }
    NOVA_REQUIRE(bad_tolerance);

    o.absolute_tolerance = 1e-4;
    o.max_attempts = 0U;
    bool bad_attempts = false;
    try { (void)solve_adaptive_heat_3d(p, o); }
    catch (const std::invalid_argument&) { bad_attempts = true; }
    NOVA_REQUIRE(bad_attempts);

    auto decaying = p;
    decaying.initial_temperature = {0.0, 0.0, 0.0, 10.0};
    decaying.diffusion.dirichlet_constraints = {{0, 0.0}, {1, 0.0}, {2, 0.0}};
    o.max_attempts = 20U;
    o.initial_step = 0.5;
    o.min_step = 0.5;
    o.max_step = 0.5;
    o.absolute_tolerance = 1e-8;
    o.relative_tolerance = 0.0;
    bool minimum_failed = false;
    try { (void)solve_adaptive_heat_3d(decaying, o); }
    catch (const std::runtime_error&) { minimum_failed = true; }
    NOVA_REQUIRE(minimum_failed);

    o.min_step = 1e-6;
    o.max_attempts = 1U;
    bool attempts_failed = false;
    try { (void)solve_adaptive_heat_3d(decaying, o); }
    catch (const std::runtime_error&) { attempts_failed = true; }
    NOVA_REQUIRE(attempts_failed);
}

NOVA_TEST(adaptive_heat_3d_observer_exception_propagates) {
    auto p = base_tetra();
    p.initial_temperature = {5.0, 5.0, 5.0, 5.0};
    AdaptiveHeat3DOptions o;
    o.on_step = [](std::size_t, double, const std::vector<double>&) {
        throw std::runtime_error("observer aborted");
    };
    bool propagated = false;
    try { (void)solve_adaptive_heat_3d(p, o); }
    catch (const std::runtime_error& e) {
        propagated = std::string(e.what()) == "observer aborted";
    }
    NOVA_REQUIRE(propagated);
}
