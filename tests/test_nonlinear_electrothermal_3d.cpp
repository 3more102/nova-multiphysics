#include "test_framework.hpp"
#include "nova/physics/nonlinear_electrothermal_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
using namespace nova::fem;
using namespace nova::physics;

bool x_face(const TetrahedralMesh3D& mesh, const BoundaryFace3D& face, double x) {
    return std::all_of(face.begin(), face.end(),
        [&](std::size_t i) { return mesh.nodes[i].x == x; });
}

NonlinearElectrothermal3DProblem plate(double left = 0.0, double right = 10.0) {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2, 2, 2);
    NonlinearElectrothermal3DProblem p;
    p.base.electrical.mesh = box.mesh;
    p.base.thermal.mesh = box.mesh;
    p.base.electrical.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    p.base.thermal.cell_diffusivity.assign(box.mesh.elements.size(), 4.0);
    p.conductivity_scale = [](double) { return 1.0; };
    for (const auto& f : box.boundary_faces) {
        if (x_face(box.mesh, f, 0.0)) {
            p.base.electrical.dirichlet_faces.push_back({
                f, [left](double, double, double) { return left; }
            });
            p.base.thermal.dirichlet_faces.push_back({
                f, [](double, double, double) { return 300.0; }
            });
        } else if (x_face(box.mesh, f, 1.0)) {
            p.base.electrical.dirichlet_faces.push_back({
                f, [right](double, double, double) { return right; }
            });
            p.base.thermal.dirichlet_faces.push_back({
                f, [](double, double, double) { return 300.0; }
            });
        }
    }
    return p;
}

template<class F> bool rejects(F&& f) {
    try { f(); } catch (const std::exception&) { return true; }
    return false;
}
} // namespace

NOVA_TEST(nonlinear_electrothermal_3d_constant_law_matches_one_way_solution) {
    const auto p = plate();
    const auto ref = solve_one_way_electrothermal_3d(p.base);
    const auto r = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE(r.iterations >= 2);
    NOVA_REQUIRE(r.max_normalized_residual <= 1.0);
    NOVA_REQUIRE_NEAR(r.fields.total_joule_power_watts, 200.0, 1e-7);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.generated_heat, 200.0, 1e-7);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.signed_balance_defect, 0.0, 2e-7);
    for (std::size_t i = 0; i < ref.thermal_solution.values.size(); ++i) {
        NOVA_REQUIRE_NEAR(r.fields.thermal_solution.values[i],
                          ref.thermal_solution.values[i], 1e-8);
    }
    NOVA_REQUIRE_NEAR(p.base.electrical.cell_diffusivity[0], 2.0, 0.0);
    NOVA_REQUIRE(p.base.thermal.cell_source.empty());
}

NOVA_TEST(nonlinear_electrothermal_3d_positive_feedback_and_voltage_gauge) {
    auto p = plate();
    p.conductivity_scale = [](double t_kelvin) {
        return 1.0 + 0.003 * (t_kelvin - 300.0);
    };
    const auto baseline = solve_one_way_electrothermal_3d(p.base);
    const auto a = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(a.converged);
    NOVA_REQUIRE(a.iterations > 2);
    NOVA_REQUIRE(a.fields.total_joule_power_watts > baseline.total_joule_power_watts);
    NOVA_REQUIRE_NEAR(a.fields.thermal_balance.signed_balance_defect, 0.0, 1e-6);
    NOVA_REQUIRE(a.fields.thermal_balance.max_free_equation_residual < 1e-6);

    auto shifted = plate(50.0, 60.0);
    shifted.conductivity_scale = p.conductivity_scale;
    const auto b = solve_nonlinear_electrothermal_3d(shifted);
    NOVA_REQUIRE(b.converged);
    NOVA_REQUIRE_NEAR(a.fields.total_joule_power_watts,
                      b.fields.total_joule_power_watts, 1e-5);
    for (std::size_t i = 0; i < a.fields.thermal_solution.values.size(); ++i) {
        NOVA_REQUIRE_NEAR(a.fields.thermal_solution.values[i],
                          b.fields.thermal_solution.values[i], 1e-6);
    }

    // The cell temperature of the converged thermal field must reconstruct
    // the material used for the last electrical solve (within Picard tolerance).
    for (std::size_t cell = 0; cell < p.base.electrical.mesh.elements.size(); ++cell) {
        double mean = 0.0;
        for (const auto node : p.base.electrical.mesh.elements[cell])
            mean += a.fields.thermal_solution.values[node] / 4.0;
        double mean_old = 0.0;
        for (const auto node : p.base.electrical.mesh.elements[cell])
            mean_old += a.last_temperature_iterate[node] / 4.0;
        NOVA_REQUIRE_NEAR(p.conductivity_scale(mean), p.conductivity_scale(mean_old), 1e-7);
    }
}

NOVA_TEST(nonlinear_electrothermal_3d_tensor_material_and_base_source) {
    auto p = plate(0.0, 0.0);
    p.base.electrical.cell_diffusivity.clear();
    p.base.electrical.cell_diffusivity_tensors.assign(
        p.base.electrical.mesh.elements.size(), {2.0, 3.0, 4.0, 0.1, 0.0, 0.0});
    p.base.thermal.cell_source.assign(p.base.thermal.mesh.elements.size(), 7.0);
    p.base.thermal.source = [](double, double, double) { return 3.0; };
    p.conductivity_scale = [](double t) { return 1.0 + 0.001 * (t - 300.0); };
    const auto r = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE_NEAR(r.fields.total_joule_power_watts, 0.0, 1e-10);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.generated_heat, 10.0, 1e-8);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.signed_balance_defect, 0.0, 1e-7);
    NOVA_REQUIRE_NEAR(p.base.electrical.cell_diffusivity_tensors.front().xx, 2.0, 0.0);
}

NOVA_TEST(nonlinear_electrothermal_3d_reports_nonconvergence_without_false_success) {
    const auto p = plate();
    NonlinearElectrothermal3DOptions options;
    options.max_iterations = 1;
    const auto first = solve_nonlinear_electrothermal_3d(p, options);
    NOVA_REQUIRE(!first.converged);
    NOVA_REQUIRE(first.iterations == 1);
    NOVA_REQUIRE(first.max_fixed_point_residual_kelvin > 1.0);
    NOVA_REQUIRE(first.max_normalized_residual > 1.0);

    // A tiny relaxation must not hide a large *raw* fixed-point defect.
    options.max_iterations = 2;
    options.relaxation = 1e-12;
    const auto slowly = solve_nonlinear_electrothermal_3d(p, options);
    NOVA_REQUIRE(!slowly.converged);
    NOVA_REQUIRE(slowly.max_fixed_point_residual_kelvin > 1.0);
}

NOVA_TEST(nonlinear_electrothermal_3d_rejects_bad_feedback_and_options) {
    auto p = plate();
    p.conductivity_scale = {};
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.conductivity_scale = [](double) { return -1.0; };
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.conductivity_scale = [](double) { return std::numeric_limits<double>::infinity(); };
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.base.electrical.diffusivity = [](double, double, double) { return 2.0; };
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.base.electrical.cell_diffusivity.pop_back();
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.initial_temperature_kelvin = -20.0;
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    auto options = NonlinearElectrothermal3DOptions{};
    options.max_iterations = 0;
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p, options); }));
    options = {};
    options.relaxation = 0.0;
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p, options); }));
    options = {};
    options.absolute_tolerance_kelvin = 0.0;
    options.relative_tolerance = 0.0;
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p, options); }));
}

NOVA_TEST(nonlinear_electrothermal_3d_variable_thermal_k_reduces_hotspot) {
    auto p = plate();
    const auto constant = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(constant.converged);
    p.thermal_conductivity_scale = [](double t_kelvin) {
        return 1.0 + 0.02 * (t_kelvin - 300.0);
    };
    const auto r = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE(r.iterations > 2);
    NOVA_REQUIRE_NEAR(r.fields.total_joule_power_watts, 200.0, 1e-7);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.generated_heat, 200.0, 1e-7);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.signed_balance_defect, 0.0, 2e-7);
    NOVA_REQUIRE(r.fields.thermal_balance.max_free_equation_residual < 2e-7);
    bool found_cooler_interior = false;
    for (std::size_t i = 0; i < p.base.thermal.mesh.nodes.size(); ++i) {
        const double x = p.base.thermal.mesh.nodes[i].x;
        if (x == 0.5 &&
            r.fields.thermal_solution.values[i] <
                constant.fields.thermal_solution.values[i] - 0.1) {
            found_cooler_interior = true;
        }
    }
    NOVA_REQUIRE(found_cooler_interior);
    NOVA_REQUIRE(r.fields.heated_thermal_problem.cell_diffusivity.size() ==
                 p.base.thermal.mesh.elements.size());
    for (std::size_t cell = 0; cell < p.base.thermal.mesh.elements.size(); ++cell) {
        double mean = 0.0;
        for (const auto vertex : p.base.thermal.mesh.elements[cell]) {
            mean += r.last_temperature_iterate[vertex] / 4.0;
        }
        NOVA_REQUIRE_NEAR(r.fields.heated_thermal_problem.cell_diffusivity[cell],
                          4.0 * p.thermal_conductivity_scale(mean), 1e-10);
    }
    NOVA_REQUIRE_NEAR(p.base.thermal.cell_diffusivity.front(), 4.0, 0.0);
}

NOVA_TEST(nonlinear_electrothermal_3d_combined_sigma_and_k_feedback) {
    auto p = plate();
    p.conductivity_scale = [](double t_kelvin) {
        return 1.0 + 0.002 * (t_kelvin - 300.0);
    };
    p.thermal_conductivity_scale = [](double t_kelvin) {
        return 1.0 + 0.01 * (t_kelvin - 300.0);
    };
    const auto r = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE(r.iterations > 2);
    NOVA_REQUIRE(r.fields.total_joule_power_watts > 200.0);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.generated_heat,
                      r.fields.total_joule_power_watts, 1e-6);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.signed_balance_defect, 0.0, 1e-6);
    NOVA_REQUIRE(r.fields.thermal_balance.max_free_equation_residual < 1e-6);

    auto shifted = plate(20.0, 30.0);
    shifted.conductivity_scale = p.conductivity_scale;
    shifted.thermal_conductivity_scale = p.thermal_conductivity_scale;
    const auto b = solve_nonlinear_electrothermal_3d(shifted);
    NOVA_REQUIRE(b.converged);
    NOVA_REQUIRE_NEAR(b.fields.total_joule_power_watts,
                      r.fields.total_joule_power_watts, 1e-5);
    for (std::size_t i = 0; i < b.fields.thermal_solution.values.size(); ++i) {
        NOVA_REQUIRE_NEAR(b.fields.thermal_solution.values[i],
                          r.fields.thermal_solution.values[i], 1e-6);
    }
}

NOVA_TEST(nonlinear_electrothermal_3d_tensor_thermal_feedback_and_additive_heat) {
    auto p = plate(0.0, 0.0);
    p.base.thermal.cell_diffusivity.clear();
    p.base.thermal.cell_diffusivity_tensors.assign(
        p.base.thermal.mesh.elements.size(), {4.0, 5.0, 6.0, 0.1, 0.0, 0.0});
    p.base.thermal.cell_source.assign(p.base.thermal.mesh.elements.size(), 7.0);
    p.base.thermal.source = [](double, double, double) { return 3.0; };
    p.thermal_conductivity_scale = [](double t_kelvin) {
        return 1.0 + 0.002 * (t_kelvin - 300.0);
    };
    const auto r = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE_NEAR(r.fields.total_joule_power_watts, 0.0, 1e-10);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.generated_heat, 10.0, 1e-8);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.signed_balance_defect, 0.0, 1e-7);
    const auto& updated = r.fields.heated_thermal_problem.cell_diffusivity_tensors;
    NOVA_REQUIRE(updated.size() == p.base.thermal.mesh.elements.size());
    for (std::size_t cell = 0; cell < updated.size(); ++cell) {
        double mean = 0.0;
        for (const auto node : p.base.thermal.mesh.elements[cell]) {
            mean += r.last_temperature_iterate[node] / 4.0;
        }
        const double multiplier = p.thermal_conductivity_scale(mean);
        NOVA_REQUIRE_NEAR(updated[cell].xx, 4.0 * multiplier, 1e-10);
        NOVA_REQUIRE_NEAR(updated[cell].yy, 5.0 * multiplier, 1e-10);
        NOVA_REQUIRE_NEAR(updated[cell].xy, 0.1 * multiplier, 1e-10);
    }
    NOVA_REQUIRE_NEAR(p.base.thermal.cell_diffusivity_tensors[0].xx, 4.0, 0.0);
}

NOVA_TEST(nonlinear_electrothermal_3d_constant_k_preserves_spatial_thermal_material) {
    auto p = plate();
    p.base.thermal.cell_diffusivity.clear();
    p.base.thermal.diffusivity = [](double, double, double) { return 4.0; };
    const auto r = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE_NEAR(r.fields.total_joule_power_watts, 200.0, 1e-7);
    NOVA_REQUIRE_NEAR(r.fields.thermal_balance.signed_balance_defect, 0.0, 2e-7);
}

NOVA_TEST(nonlinear_electrothermal_3d_rejects_invalid_thermal_feedback) {
    auto p = plate();
    p.thermal_conductivity_scale = [](double) { return 0.0; };
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.thermal_conductivity_scale = [](double) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.thermal_conductivity_scale = [](double) { return -1.0; };
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.thermal_conductivity_scale = [](double) { return 1.0; };
    p.base.thermal.cell_diffusivity.pop_back();
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.thermal_conductivity_scale = [](double) { return 1.0; };
    p.base.thermal.diffusivity = [](double, double, double) { return 4.0; };
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
    p = plate();
    p.thermal_conductivity_scale = [](double) { return 1.0; };
    p.base.thermal.cell_diffusivity_tensors.assign(
        p.base.thermal.mesh.elements.size(), {4.0, 4.0, 4.0, 0.0, 0.0, 0.0});
    NOVA_REQUIRE(rejects([&] { (void)solve_nonlinear_electrothermal_3d(p); }));
}

NOVA_TEST(nonlinear_electrothermal_3d_history_records_every_solved_step) {
    auto p = plate();
    const auto r = solve_nonlinear_electrothermal_3d(p);
    NOVA_REQUIRE(r.converged);
    NOVA_REQUIRE(r.history.size() == r.iterations);
    NOVA_REQUIRE(r.history.size() >= 2);
    for (std::size_t i = 0; i < r.history.size(); ++i) {
        const auto& h = r.history[i];
        NOVA_REQUIRE(h.iteration == i + 1);
        NOVA_REQUIRE_NEAR(h.total_joule_power_watts, 200.0, 1e-7);
        NOVA_REQUIRE_NEAR(h.signed_thermal_balance_defect_watts, 0.0, 2e-7);
        NOVA_REQUIRE(h.min_solved_temperature_kelvin >= 300.0 - 1e-8);
        NOVA_REQUIRE(h.max_solved_temperature_kelvin >= h.min_solved_temperature_kelvin);
        NOVA_REQUIRE(h.max_normalized_residual >= 0.0);
        NOVA_REQUIRE(h.max_fixed_point_residual_kelvin >= 0.0);
        if (i + 1 < r.history.size()) {
            NOVA_REQUIRE_NEAR(h.relaxation_applied, 0.7, 0.0);
        } else {
            NOVA_REQUIRE_NEAR(h.relaxation_applied, 0.0, 0.0);
        }
    }
    NOVA_REQUIRE_NEAR(r.history.back().max_fixed_point_residual_kelvin,
                      r.max_fixed_point_residual_kelvin, 0.0);
    NOVA_REQUIRE_NEAR(r.history.back().max_normalized_residual,
                      r.max_normalized_residual, 0.0);
}

NOVA_TEST(nonlinear_electrothermal_3d_adaptive_damping_is_bounded_and_traceable) {
    auto p = plate();
    // Strong but positive temperature coefficient; restrict the iteration
    // horizon so this does not assume global nonlinear convergence.
    p.conductivity_scale = [](double t) {
        return 1.0 + 0.8 * std::max(0.0, t - 300.0);
    };
    NonlinearElectrothermal3DOptions options;
    options.adaptive_relaxation = true;
    options.relaxation = 1.0;
    options.minimum_relaxation = 0.125;
    options.max_iterations = 5;
    const auto r = solve_nonlinear_electrothermal_3d(p, options);
    NOVA_REQUIRE(r.history.size() == r.iterations);
    NOVA_REQUIRE(!r.history.empty());
    NOVA_REQUIRE(r.iterations <= options.max_iterations);
    NOVA_REQUIRE_NEAR(r.history.back().relaxation_applied, 0.0, 0.0);
    for (std::size_t i = 0; i < r.history.size(); ++i) {
        const auto& h = r.history[i];
        NOVA_REQUIRE(h.iteration == i + 1);
        NOVA_REQUIRE(std::isfinite(h.max_normalized_residual));
        NOVA_REQUIRE(std::isfinite(h.max_fixed_point_residual_kelvin));
        NOVA_REQUIRE(std::isfinite(h.total_joule_power_watts));
        NOVA_REQUIRE(std::isfinite(h.signed_thermal_balance_defect_watts));
        if (i + 1 == r.history.size()) continue;
        NOVA_REQUIRE(h.relaxation_applied >= options.minimum_relaxation);
        NOVA_REQUIRE(h.relaxation_applied <= options.relaxation);
        if (i == 0) {
            NOVA_REQUIRE_NEAR(h.relaxation_applied, options.relaxation, 0.0);
        } else {
            const auto& prev = r.history[i - 1];
            const long double curr =
                static_cast<long double>(h.max_normalized_residual);
            const long double prior =
                static_cast<long double>(prev.max_normalized_residual);
            if (curr > 1.05L * prior) {
                NOVA_REQUIRE_NEAR(
                    h.relaxation_applied,
                    std::max(options.minimum_relaxation, 0.5 * prev.relaxation_applied),
                    1e-14);
            } else if (curr < 0.8L * prior) {
                NOVA_REQUIRE_NEAR(
                    h.relaxation_applied,
                    std::min(options.relaxation, 1.2 * prev.relaxation_applied),
                    1e-14);
            } else {
                NOVA_REQUIRE_NEAR(h.relaxation_applied, prev.relaxation_applied, 0.0);
            }
        }
    }
}

NOVA_TEST(nonlinear_electrothermal_3d_history_on_iteration_cap_and_invalid_damping) {
    const auto p = plate();
    NonlinearElectrothermal3DOptions options;
    options.max_iterations = 1;
    options.adaptive_relaxation = true;
    const auto r = solve_nonlinear_electrothermal_3d(p, options);
    NOVA_REQUIRE(!r.converged);
    NOVA_REQUIRE(r.iterations == 1);
    NOVA_REQUIRE(r.history.size() == 1);
    NOVA_REQUIRE_NEAR(r.history.front().relaxation_applied, 0.0, 0.0);
    NOVA_REQUIRE_NEAR(r.history.front().max_normalized_residual,
                      r.max_normalized_residual, 0.0);
    options.minimum_relaxation = 0.8;
    NOVA_REQUIRE(rejects([&] {
        (void)solve_nonlinear_electrothermal_3d(p, options);
    }));
    options.minimum_relaxation = 0.0;
    NOVA_REQUIRE(rejects([&] {
        (void)solve_nonlinear_electrothermal_3d(p, options);
    }));
    options.minimum_relaxation = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] {
        (void)solve_nonlinear_electrothermal_3d(p, options);
    }));
}
