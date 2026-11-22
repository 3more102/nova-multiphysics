#include "test_framework.hpp"
#include "nova/physics/electrothermal_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {
using namespace nova::fem;
using namespace nova::physics;

bool at_x(const TetrahedralMesh3D& mesh, const BoundaryFace3D& face, double x) {
    return std::all_of(face.begin(), face.end(),
        [&](std::size_t i) { return mesh.nodes[i].x == x; });
}

Electrothermal3DProblem make_plate(double left_v, double right_v) {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2, 2, 2);
    Electrothermal3DProblem p;
    p.electrical.mesh = box.mesh;
    p.thermal.mesh = box.mesh;
    p.electrical.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    p.thermal.cell_diffusivity.assign(box.mesh.elements.size(), 4.0);
    for (const auto& f : box.boundary_faces) {
        if (at_x(box.mesh, f, 0.0)) {
            p.electrical.dirichlet_faces.push_back({
                f, [left_v](double, double, double) { return left_v; }
            });
            p.thermal.dirichlet_faces.push_back({
                f, [](double, double, double) { return 300.0; }
            });
        } else if (at_x(box.mesh, f, 1.0)) {
            p.electrical.dirichlet_faces.push_back({
                f, [right_v](double, double, double) { return right_v; }
            });
            p.thermal.dirichlet_faces.push_back({
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

NOVA_TEST(electrothermal_3d_dc_parallel_plate_and_thermal_power_conservation) {
    const auto p = make_plate(0.0, 10.0);
    const auto r = solve_one_way_electrothermal_3d(p);
    NOVA_REQUIRE(r.electrical_solution.linear_solver.converged);
    NOVA_REQUIRE(r.thermal_solution.linear_solver.converged);
    NOVA_REQUIRE(r.joule_density.size() == p.electrical.mesh.elements.size());
    NOVA_REQUIRE_NEAR(r.total_joule_power_watts, 200.0, 1e-7);
    NOVA_REQUIRE_NEAR(r.thermal_balance.generated_heat, 200.0, 1e-7);
    NOVA_REQUIRE_NEAR(r.thermal_balance.signed_balance_defect, 0.0, 2e-7);
    NOVA_REQUIRE(r.thermal_balance.max_free_equation_residual < 2e-7);
    NOVA_REQUIRE_NEAR(r.thermal_balance.constraint_outward_heat, 200.0, 2e-7);
    for (std::size_t i = 0; i < r.joule_density.size(); ++i) {
        NOVA_REQUIRE_NEAR(r.electric_field[i].x, -10.0, 1e-8);
        NOVA_REQUIRE_NEAR(r.current_density[i].x, -20.0, 1e-8);
        NOVA_REQUIRE_NEAR(r.current_density[i].y, 0.0, 1e-8);
        NOVA_REQUIRE_NEAR(r.joule_density[i], 200.0, 1e-7);
    }
    for (std::size_t i = 0; i < p.electrical.mesh.nodes.size(); ++i) {
        const double x = p.electrical.mesh.nodes[i].x;
        NOVA_REQUIRE_NEAR(r.electrical_solution.values[i], 10.0*x, 1e-8);
        if (x == 0.0 || x == 1.0) {
            NOVA_REQUIRE_NEAR(r.thermal_solution.values[i], 300.0, 1e-8);
        }
        if (x == 0.5) {
            NOVA_REQUIRE_NEAR(r.thermal_solution.values[i], 306.25, 0.5);
        }
    }
    NOVA_REQUIRE(p.thermal.cell_source.empty());
}

NOVA_TEST(electrothermal_3d_spd_tensor_gauge_invariance_and_joule_density) {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 1, 1, 1);
    auto make = [&](double shift) {
        Electrothermal3DProblem p;
        p.electrical.mesh = box.mesh;
        p.thermal.mesh = box.mesh;
        p.electrical.cell_diffusivity_tensors.assign(
            box.mesh.elements.size(), {2, 3, 4, 0.25, 0, 0});
        p.thermal.cell_diffusivity.assign(box.mesh.elements.size(), 3.0);
        for (const auto& face : box.boundary_faces) {
            p.electrical.dirichlet_faces.push_back({
                face, [shift](double x, double, double) { return x + shift; }
            });
            p.thermal.dirichlet_faces.push_back({
                face, [](double, double, double) { return 300.0; }
            });
        }
        return p;
    };
    const auto a = solve_one_way_electrothermal_3d(make(0));
    const auto b = solve_one_way_electrothermal_3d(make(4));
    NOVA_REQUIRE_NEAR(a.total_joule_power_watts, 2.0, 1e-10);
    NOVA_REQUIRE_NEAR(b.total_joule_power_watts, 2.0, 1e-10);
    NOVA_REQUIRE_NEAR(a.thermal_balance.generated_heat, 2.0, 1e-10);
    NOVA_REQUIRE_NEAR(a.thermal_balance.constraint_outward_heat, 2.0, 1e-9);
    for (std::size_t i = 0; i < a.joule_density.size(); ++i) {
        NOVA_REQUIRE_NEAR(a.joule_density[i], 2.0, 1e-10);
        NOVA_REQUIRE_NEAR(a.current_density[i].x, -2.0, 1e-10);
        NOVA_REQUIRE_NEAR(a.current_density[i].y, -0.25, 1e-10);
        NOVA_REQUIRE_NEAR(a.joule_density[i], b.joule_density[i], 1e-10);
    }
}

NOVA_TEST(electrothermal_3d_additive_base_heat_sources_and_zero_voltage_gradient) {
    auto p = make_plate(5.0, 5.0);
    p.thermal.cell_source.assign(p.thermal.mesh.elements.size(), 7.0);
    p.thermal.source = [](double, double, double) { return 3.0; };
    const auto r = solve_one_way_electrothermal_3d(p);
    NOVA_REQUIRE_NEAR(r.total_joule_power_watts, 0.0, 1e-11);
    NOVA_REQUIRE_NEAR(r.thermal_balance.generated_heat, 10.0, 1e-9);
    NOVA_REQUIRE_NEAR(r.thermal_balance.constraint_outward_heat, 10.0, 1e-7);
    for (double q : r.joule_density) NOVA_REQUIRE_NEAR(q, 0.0, 1e-11);
    for (double q : r.heated_thermal_problem.cell_source)
        NOVA_REQUIRE_NEAR(q, 7.0, 1e-11);
    for (double q : p.thermal.cell_source) NOVA_REQUIRE_NEAR(q, 7.0, 1e-11);
}

NOVA_TEST(electrothermal_3d_rejects_misaligned_mesh_and_invalid_materials) {
    auto p = make_plate(0.0, 1.0);
    p.thermal.mesh.nodes.front().x += 0.001;
    NOVA_REQUIRE(rejects([&] { (void)solve_one_way_electrothermal_3d(p); }));
    p = make_plate(0.0, 1.0);
    std::swap(p.thermal.mesh.elements[0], p.thermal.mesh.elements[1]);
    NOVA_REQUIRE(rejects([&] { (void)solve_one_way_electrothermal_3d(p); }));
    p = make_plate(0.0, 1.0);
    p.electrical.cell_diffusivity.front() = -1.0;
    NOVA_REQUIRE(rejects([&] { (void)solve_one_way_electrothermal_3d(p); }));
    p = make_plate(0.0, 1.0);
    p.thermal.cell_source = {1.0};
    NOVA_REQUIRE(rejects([&] { (void)solve_one_way_electrothermal_3d(p); }));
    p = make_plate(0.0, 1.0);
    p.thermal.cell_source.assign(
        p.thermal.mesh.elements.size(), std::numeric_limits<double>::infinity());
    NOVA_REQUIRE(rejects([&] { (void)solve_one_way_electrothermal_3d(p); }));
}
