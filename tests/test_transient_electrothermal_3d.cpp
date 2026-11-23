#include "test_framework.hpp"
#include "nova/physics/transient_electrothermal_3d.hpp"
#include "nova/physics/adaptive_transient_electrothermal_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

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
    const auto r = solve_adaptive_transient_electrothermal_3d(p, o);
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
