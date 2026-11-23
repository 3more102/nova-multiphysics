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
