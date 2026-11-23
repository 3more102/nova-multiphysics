#include "test_framework.hpp"
#include "nova/physics/transient_electrothermal_3d.hpp"

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
