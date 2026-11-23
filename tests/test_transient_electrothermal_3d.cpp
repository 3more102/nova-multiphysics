#include "test_framework.hpp"
#include "nova/physics/transient_electrothermal_3d.hpp"
#include "nova/physics/adaptive_transient_electrothermal_3d.hpp"
#include "nova/io/thermal_checkpoint_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>

namespace {
using namespace nova::fem;
using namespace nova::physics;

bool face_at_x(const TetrahedralMesh3D& mesh,
               const BoundaryFace3D& f, double x) {
    return std::all_of(f.begin(), f.end(),
        [&](std::size_t i) { return mesh.nodes[i].x == x; });
}

TransientElectrothermal3DProblem make_unit_coupling() {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 1, 1, 1);
    TransientElectrothermal3DProblem p;
    p.electrical.mesh = box.mesh;
    p.electrical.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    p.thermal.diffusion.mesh = box.mesh;
    p.thermal.diffusion.cell_diffusivity.assign(box.mesh.elements.size(), 4.0);
    p.thermal.cell_volumetric_capacity.assign(box.mesh.elements.size(), 2.0);
    p.thermal.initial_temperature.assign(box.mesh.nodes.size(), 300.0);
    for (const auto& f : box.boundary_faces) {
        if (face_at_x(box.mesh, f, 0.0) ||
            face_at_x(box.mesh, f, 1.0)) {
            p.electrical.dirichlet_faces.push_back({
                f, [](double x, double, double) { return 10.0 * x; }
            });
        }
    }
    return p;
}

template <typename F> bool rejects(F&& f) {
    try { f(); } catch (const std::exception&) { return true; }
    return false;
}
} // namespace

NOVA_TEST(transient_electrothermal_3d_uniform_joule_ramp_and_energy_budget) {
    const auto p = make_unit_coupling();
    TransientElectrothermal3DOptions options;
    options.thermal_time.time_step = 0.02;
    options.thermal_time.steps = 4;
    const auto r = solve_transient_electrothermal_3d(p, options);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE(r.failed_step == 0);
    NOVA_REQUIRE(r.accepted_steps.size() == 4);
    NOVA_REQUIRE_NEAR(r.final_time, 0.08, 1e-14);
    NOVA_REQUIRE_NEAR(r.last_electrical.total_joule_power_watts, 200.0, 1e-7);
    for (const double t : r.temperature) {
        NOVA_REQUIRE_NEAR(t, 308.0, 2e-8);
    }
    for (std::size_t j = 0; j < r.accepted_steps.size(); ++j) {
        const auto& s = r.accepted_steps[j];
        NOVA_REQUIRE(s.coupling_iterations >= 2);
        NOVA_REQUIRE(s.max_normalized_residual <= 1.0);
        NOVA_REQUIRE_NEAR(s.time, (j + 1) * 0.02, 1e-14);
        NOVA_REQUIRE_NEAR(s.joule_power_watts, 200.0, 1e-7);
        NOVA_REQUIRE_NEAR(s.heat_budget.source_energy, 4.0, 1e-8);
        NOVA_REQUIRE_NEAR(s.heat_budget.heat_content,
                          2.0 * (302.0 + 2.0 * j), 1e-7);
        NOVA_REQUIRE_NEAR(s.heat_budget.balance_residual, 0.0, 1e-8);
    }
    for (double v : p.thermal.initial_temperature)
        NOVA_REQUIRE_NEAR(v, 300.0, 0.0);
    NOVA_REQUIRE(p.thermal.diffusion.cell_source.empty());
}

NOVA_TEST(transient_electrothermal_3d_nonlin_sigma_and_k_match_implicit_uniform_ode) {
    auto p = make_unit_coupling();
    p.electrical_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };
    p.thermal_conductivity_scale = [](double t) {
        return 1.0 + 0.02 * (t - 300.0);
    };
    TransientElectrothermal3DOptions options;
    options.thermal_time.time_step = 0.05;
    options.thermal_time.steps = 3;
    const auto r = solve_transient_electrothermal_3d(p, options);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE(r.accepted_steps.size() == 3);
    double exact_delta = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        // C=2 J/(m3 K), P(T)=200*(1+0.01*(T-300)) W/m3.
        // Backward Euler: delta_next=(delta_old+100*dt)/(1-dt).
        exact_delta = (exact_delta + 5.0) / 0.95;
        const auto& s = r.accepted_steps[i];
        const double exact_power = 200.0 * (1.0 + 0.01 * exact_delta);
        NOVA_REQUIRE_NEAR(s.joule_power_watts, exact_power, 1e-5);
        NOVA_REQUIRE_NEAR(s.heat_budget.source_energy,
                          options.thermal_time.time_step * exact_power, 1e-6);
        NOVA_REQUIRE_NEAR(s.heat_budget.balance_residual, 0.0, 1e-7);
        NOVA_REQUIRE(s.max_normalized_residual <= 1.0);
    }
    for (const double t : r.temperature) {
        NOVA_REQUIRE_NEAR(t, 300.0 + exact_delta, 1e-6);
    }
    NOVA_REQUIRE_NEAR(p.electrical.cell_diffusivity.front(), 2.0, 0.0);
    NOVA_REQUIRE_NEAR(p.thermal.diffusion.cell_diffusivity.front(), 4.0, 0.0);
}

NOVA_TEST(transient_electrothermal_3d_tensor_and_additive_heat) {
    auto p = make_unit_coupling();
    p.electrical.cell_diffusivity.clear();
    p.electrical.cell_diffusivity_tensors.assign(
        p.electrical.mesh.elements.size(), {2, 3, 4, 0, 0, 0});
    p.thermal.diffusion.cell_diffusivity.clear();
    p.thermal.diffusion.cell_diffusivity_tensors.assign(
        p.thermal.diffusion.mesh.elements.size(), {4, 5, 6, 0, 0, 0});
    p.thermal.diffusion.cell_source.assign(p.thermal.diffusion.mesh.elements.size(), 10.0);
    p.electrical_conductivity_scale = [](double) { return 1.0; };
    p.thermal_conductivity_scale = [](double) { return 1.0; };
    TransientElectrothermal3DOptions options;
    options.thermal_time.time_step = 0.01;
    options.thermal_time.steps = 2;
    const auto r = solve_transient_electrothermal_3d(p, options);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE_NEAR(r.last_electrical.total_joule_power_watts, 200.0, 1e-7);
    for (const auto t : r.temperature)
        NOVA_REQUIRE_NEAR(t, 302.1, 1e-7);
    for (const auto& s : r.accepted_steps) {
        NOVA_REQUIRE_NEAR(s.heat_budget.source_energy, 2.1, 1e-7);
        NOVA_REQUIRE_NEAR(s.heat_budget.balance_residual, 0.0, 1e-7);
    }
}

NOVA_TEST(transient_electrothermal_3d_rejects_nonconverged_step_without_advancing) {
    auto p = make_unit_coupling();
    TransientElectrothermal3DOptions options;
    options.max_coupling_iterations = 1;
    options.thermal_time.time_step = 0.05;
    options.thermal_time.steps = 2;
    const auto r = solve_transient_electrothermal_3d(p, options);
    NOVA_REQUIRE(!r.converged);
    NOVA_REQUIRE(r.failed_step == 1);
    NOVA_REQUIRE(r.accepted_steps.empty());
    NOVA_REQUIRE(r.failed_max_normalized_residual > 1.0);
    NOVA_REQUIRE_NEAR(r.final_time, 0.0, 0.0);
    for (const auto t : r.temperature)
        NOVA_REQUIRE_NEAR(t, 300.0, 0.0);
}

NOVA_TEST(transient_electrothermal_3d_rejects_incompatible_materials_and_options) {
    auto p = make_unit_coupling();
    p.thermal.diffusion.mesh.nodes.front().x += 0.01;
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_electrothermal_3d(p); }));
    p = make_unit_coupling();
    p.thermal.diffusion.cell_source = {1.0};
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_electrothermal_3d(p); }));
    p = make_unit_coupling();
    p.electrical_conductivity_scale = [](double) { return -1.0; };
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_electrothermal_3d(p); }));
    p = make_unit_coupling();
    p.thermal_conductivity_scale = [](double) {
        return std::numeric_limits<double>::infinity();
    };
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_electrothermal_3d(p); }));
    p = make_unit_coupling();
    p.thermal_conductivity_scale = [](double) { return 1.0; };
    p.thermal.diffusion.diffusivity = [](double, double, double) { return 4.0; };
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_electrothermal_3d(p); }));
    p = make_unit_coupling();
    auto opts = TransientElectrothermal3DOptions{};
    opts.thermal_time.scheme = TransientHeat3DScheme::CrankNicolson;
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_electrothermal_3d(p, opts); }));
    opts = {};
    opts.thermal_time.on_step = [](std::size_t, double, const std::vector<double>&) {};
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_electrothermal_3d(p, opts); }));
    opts = {};
    opts.max_coupling_iterations = 0;
    NOVA_REQUIRE(rejects([&] { (void)solve_transient_electrothermal_3d(p, opts); }));
}

NOVA_TEST(adaptive_transient_electrothermal_3d_constant_joule_ramp_energy) {
    const auto p = make_unit_coupling();
    AdaptiveTransientElectrothermal3DOptions o;
    o.initial_time_step = 0.10;
    o.maximum_time_step = 0.10;
    o.end_time = 0.14;
    const auto r = solve_adaptive_transient_electrothermal_3d(p, o);
    NOVA_REQUIRE(r.completed);
    NOVA_REQUIRE(r.stop_reason ==
                 AdaptiveTransientElectrothermal3DStopReason::Completed);
    NOVA_REQUIRE(r.rejected_steps == 0);
    NOVA_REQUIRE(r.accepted_steps.size() == 2);
    NOVA_REQUIRE(r.attempts.size() == 2);
    NOVA_REQUIRE_NEAR(r.final_time, 0.14, 1e-14);
    for (const double t : r.temperature) {
        NOVA_REQUIRE_NEAR(t, 314.0, 2e-8);
    }
    double integrated_joules = 0.0;
    double elapsed = 0.0;
    for (const auto& step : r.accepted_steps) {
        NOVA_REQUIRE(step.max_normalized_local_error <= 1.0);
        NOVA_REQUIRE(step.nonlinear_iterations >= 3U);
        NOVA_REQUIRE_NEAR(step.joule_power_watts, 200.0, 1e-7);
        NOVA_REQUIRE_NEAR(step.heat_budget.source_energy,
                          200.0 * step.time_step, 2e-8);
        NOVA_REQUIRE_NEAR(step.heat_budget.balance_residual, 0.0, 1e-8);
        elapsed += step.time_step;
        integrated_joules += step.heat_budget.source_energy;
    }
    NOVA_REQUIRE_NEAR(elapsed, 0.14, 1e-14);
    NOVA_REQUIRE_NEAR(integrated_joules, 28.0, 2e-8);
    NOVA_REQUIRE_NEAR(r.accepted_steps.back().heat_budget.heat_content,
                      2.0 * 314.0, 2e-7);
    for (const auto t : p.thermal.initial_temperature)
        NOVA_REQUIRE_NEAR(t, 300.0, 0.0);
}

NOVA_TEST(adaptive_transient_electrothermal_3d_rejects_and_improves_nonlinear_ode) {
    auto p = make_unit_coupling();
    p.electrical_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };
    p.thermal_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };
    AdaptiveTransientElectrothermal3DOptions o;
    o.end_time = 0.3;
    o.initial_time_step = 0.3;
    o.maximum_time_step = 0.3;
    o.minimum_time_step = 1e-5;
    o.absolute_error_kelvin = 0.05;
    o.relative_error = 0.0;
    const auto r = solve_adaptive_transient_electrothermal_3d(p, o);
    NOVA_REQUIRE(r.completed);
    NOVA_REQUIRE(r.rejected_steps >= 1U);
    NOVA_REQUIRE(r.accepted_steps.size() >= 3U);
    NOVA_REQUIRE(r.attempts.size() == r.accepted_steps.size() + r.rejected_steps);
    NOVA_REQUIRE(r.attempts.front().nonlinear_converged);
    NOVA_REQUIRE(!r.attempts.front().accepted);
    NOVA_REQUIRE(r.attempts.front().max_normalized_local_error > 1.0);
    NOVA_REQUIRE_NEAR(r.final_time, 0.3, 1e-12);
    const double exact_t = 300.0 + 100.0 * (std::exp(0.3) - 1.0);
    for (const auto t : r.temperature) {
        NOVA_REQUIRE(std::abs(t - exact_t) < 1.5);
        // Compare with one coarse BE step: D=100*0.3/(1-0.3).
        NOVA_REQUIRE(std::abs(t - exact_t) <
                     std::abs(300.0 + 30.0 / 0.7 - exact_t));
    }
    for (const auto& step : r.accepted_steps) {
        NOVA_REQUIRE(step.max_normalized_local_error <= 1.0);
        NOVA_REQUIRE(step.heat_budget.source_energy > 0.0);
        NOVA_REQUIRE_NEAR(step.heat_budget.balance_residual, 0.0, 1e-7);
    }
}

NOVA_TEST(adaptive_transient_electrothermal_3d_nonlinear_failure_rolls_back) {
    const auto p = make_unit_coupling();
    AdaptiveTransientElectrothermal3DOptions o;
    o.coupling.max_coupling_iterations = 1;
    o.initial_time_step = 0.02;
    o.minimum_time_step = 0.005;
    o.end_time = 0.10;
    std::size_t observer_calls = 0U;
    o.on_accepted_step = [&](const auto&, const auto&) { ++observer_calls; };
    const auto r = solve_adaptive_transient_electrothermal_3d(p, o);
    NOVA_REQUIRE(observer_calls == 0U);
    NOVA_REQUIRE(!r.completed);
    NOVA_REQUIRE(r.stop_reason ==
                 AdaptiveTransientElectrothermal3DStopReason::MinimumTimeStep);
    NOVA_REQUIRE(r.accepted_steps.empty());
    NOVA_REQUIRE(r.rejected_steps >= 3U);
    NOVA_REQUIRE(r.attempts.size() == r.rejected_steps);
    NOVA_REQUIRE_NEAR(r.final_time, 0.0, 0.0);
    for (const auto& attempt : r.attempts) {
        NOVA_REQUIRE(!attempt.accepted);
        NOVA_REQUIRE(!attempt.nonlinear_converged);
        NOVA_REQUIRE(std::isinf(attempt.max_normalized_local_error));
    }
    for (const auto t : r.temperature)
        NOVA_REQUIRE_NEAR(t, 300.0, 0.0);
}

NOVA_TEST(adaptive_transient_electrothermal_3d_attempt_and_step_caps) {
    auto p = make_unit_coupling();
    p.electrical_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };
    AdaptiveTransientElectrothermal3DOptions o;
    o.initial_time_step = 0.3;
    o.maximum_time_step = 0.3;
    o.end_time = 0.3;
    o.absolute_error_kelvin = 0.05;
    o.relative_error = 0.0;
    o.maximum_attempts = 1U;
    const auto capped = solve_adaptive_transient_electrothermal_3d(p, o);
    NOVA_REQUIRE(!capped.completed);
    NOVA_REQUIRE(capped.stop_reason ==
                 AdaptiveTransientElectrothermal3DStopReason::AttemptLimit);
    NOVA_REQUIRE(capped.rejected_steps == 1U);
    NOVA_REQUIRE(capped.accepted_steps.empty());
    NOVA_REQUIRE_NEAR(capped.final_time, 0.0, 0.0);
    p = make_unit_coupling();
    o = {};
    o.initial_time_step = 0.1;
    o.maximum_time_step = 0.1;
    o.end_time = 0.15;
    o.maximum_accepted_steps = 1U;
    const auto capped_steps = solve_adaptive_transient_electrothermal_3d(p, o);
    NOVA_REQUIRE(!capped_steps.completed);
    NOVA_REQUIRE(capped_steps.stop_reason ==
                 AdaptiveTransientElectrothermal3DStopReason::AcceptedStepLimit);
    NOVA_REQUIRE(capped_steps.accepted_steps.size() == 1U);
    NOVA_REQUIRE_NEAR(capped_steps.final_time, 0.1, 1e-14);
    for (const auto t : capped_steps.temperature)
        NOVA_REQUIRE_NEAR(t, 310.0, 2e-8);
}

NOVA_TEST(adaptive_transient_electrothermal_3d_final_short_step_and_options) {
    const auto p = make_unit_coupling();
    AdaptiveTransientElectrothermal3DOptions o;
    o.initial_time_step = 0.05;
    o.maximum_time_step = 0.05;
    o.minimum_time_step = 0.01;
    o.end_time = 0.002;  // final physical interval is below minimum step
    const auto r = solve_adaptive_transient_electrothermal_3d(p, o);
    NOVA_REQUIRE(r.completed);
    NOVA_REQUIRE(r.accepted_steps.size() == 1U);
    NOVA_REQUIRE_NEAR(r.final_time, 0.002, 1e-14);
    for (const auto t : r.temperature)
        NOVA_REQUIRE_NEAR(t, 300.2, 2e-8);
    o = {};
    o.minimum_time_step = 0.2;
    NOVA_REQUIRE(rejects([&] {
        (void)solve_adaptive_transient_electrothermal_3d(p, o);
    }));
    o = {};
    o.absolute_error_kelvin = 0.0;
    o.relative_error = 0.0;
    NOVA_REQUIRE(rejects([&] {
        (void)solve_adaptive_transient_electrothermal_3d(p, o);
    }));
    o = {};
    o.end_time = 0.0;
    NOVA_REQUIRE(rejects([&] {
        (void)solve_adaptive_transient_electrothermal_3d(p, o);
    }));
    o = {};
    o.maximum_attempts = 0U;
    NOVA_REQUIRE(rejects([&] {
        (void)solve_adaptive_transient_electrothermal_3d(p, o);
    }));
}

NOVA_TEST(adaptive_transient_electrothermal_checkpoint_observer_and_restart) {
    // The nonlinear manufactured ODE starts with a rejected large step.
    // Only committed two-half-step solutions may be observed/checkpointed.
    auto p = make_unit_coupling();
    p.electrical_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };
    p.thermal_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };

    AdaptiveTransientElectrothermal3DOptions options;
    options.end_time = 0.3;
    options.initial_time_step = 0.3;
    options.maximum_time_step = 0.3;
    options.minimum_time_step = 1.0e-5;
    options.absolute_error_kelvin = 0.05;
    options.relative_error = 0.0;

    std::ostringstream checkpoint_output;
    double saved_next_step = 0.0;
    std::vector<double> observed_times;
    options.on_accepted_step = [&](const AdaptiveTransientElectrothermal3DStep& step,
                                   const std::vector<double>& temperature) {
        NOVA_REQUIRE(temperature.size() == p.electrical.mesh.nodes.size());
        NOVA_REQUIRE(step.time > 0.0);
        NOVA_REQUIRE(step.next_time_step >= options.minimum_time_step);
        NOVA_REQUIRE(step.next_time_step <= options.maximum_time_step);
        NOVA_REQUIRE(step.max_normalized_local_error <= 1.0);
        observed_times.push_back(step.time);
        if (observed_times.size() == 1U) {
            nova::io::write_transient_heat_checkpoint_3d(
                checkpoint_output, p.electrical.mesh, step.time, temperature);
            saved_next_step = step.next_time_step;
        }
    };

    const auto continuous = solve_adaptive_transient_electrothermal_3d(p, options);
    NOVA_REQUIRE(continuous.completed);
    NOVA_REQUIRE(continuous.rejected_steps >= 1U);
    NOVA_REQUIRE(continuous.accepted_steps.size() >= 3U);
    NOVA_REQUIRE(observed_times.size() == continuous.accepted_steps.size());
    NOVA_REQUIRE(continuous.attempts.size() ==
                 continuous.accepted_steps.size() + continuous.rejected_steps);
    for (std::size_t i = 0; i < observed_times.size(); ++i) {
        NOVA_REQUIRE_NEAR(observed_times[i], continuous.accepted_steps[i].time, 0.0);
    }

    std::istringstream checkpoint_input(checkpoint_output.str());
    const auto saved = nova::io::read_transient_heat_checkpoint_3d(
        checkpoint_input, p.electrical.mesh);
    NOVA_REQUIRE_NEAR(saved.time, observed_times.front(), 0.0);
    NOVA_REQUIRE(saved_next_step > 0.0);
    auto restarted_problem = p;
    restarted_problem.thermal.initial_temperature = saved.temperature;
    auto restarted_options = options;
    restarted_options.on_accepted_step = {};
    restarted_options.coupling.thermal_time.start_time = saved.time;
    restarted_options.initial_time_step = saved_next_step;
    const auto restarted = solve_adaptive_transient_electrothermal_3d(
        restarted_problem, restarted_options);
    NOVA_REQUIRE(restarted.completed);
    NOVA_REQUIRE_NEAR(restarted.final_time, continuous.final_time, 1e-14);
    NOVA_REQUIRE(restarted.accepted_steps.size() + 1U ==
                 continuous.accepted_steps.size());
    for (std::size_t i = 0; i < restarted.accepted_steps.size(); ++i) {
        NOVA_REQUIRE_NEAR(restarted.accepted_steps[i].time,
                          continuous.accepted_steps[i + 1U].time, 1e-12);
        NOVA_REQUIRE_NEAR(restarted.accepted_steps[i].next_time_step,
                          continuous.accepted_steps[i + 1U].next_time_step, 1e-12);
    }
    for (std::size_t i = 0; i < restarted.temperature.size(); ++i) {
        NOVA_REQUIRE_NEAR(restarted.temperature[i], continuous.temperature[i], 1e-8);
    }
    // Original problem and checkpoint are separate, immutable restart inputs.
    for (double t : p.thermal.initial_temperature) {
        NOVA_REQUIRE_NEAR(t, 300.0, 0.0);
    }
}
