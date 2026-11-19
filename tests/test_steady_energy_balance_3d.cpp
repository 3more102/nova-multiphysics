#include "test_framework.hpp"

#include "nova/fem/steady_energy_balance_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace nova::fem;

bool face_at_x(const TetrahedralMesh3D& mesh, const BoundaryFace3D& face, double x) {
    for (const auto node : face) {
        if (mesh.nodes[node].x != x) return false;
    }
    return true;
}

ScalarDiffusion3DProblem unit_tetrahedron() {
    ScalarDiffusion3DProblem p;
    p.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    p.diffusivity = [](double, double, double) { return 2.0; };
    p.source = [](double, double, double) { return 0.0; };
    return p;
}

template<class F>
bool rejects(F&& operation) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}
} // namespace

NOVA_TEST(energy_balance_recovers_heat_generation_from_dirichlet_reactions) {
    auto p = unit_tetrahedron();
    p.source = [](double, double, double) { return 20.0; };
    p.cell_source = {40.0};
    for (std::size_t i = 0; i < 4; ++i) p.dirichlet_constraints.push_back({i, 300.0});
    const auto solution = solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto budget = evaluate_steady_energy_balance_3d(p, solution.values);
    // Volume = 1/6 m^3; (20+40) W/m^3 = 10 W, even though raw P1 flux = 0.
    NOVA_REQUIRE_NEAR(budget.volumetric_generation, 10.0, 1e-12);
    NOVA_REQUIRE_NEAR(budget.dirichlet_outward, 10.0, 1e-10);
    NOVA_REQUIRE_NEAR(budget.neumann_outward, 0.0, 1e-14);
    NOVA_REQUIRE_NEAR(budget.robin_outward, 0.0, 1e-14);
    NOVA_REQUIRE_NEAR(budget.imbalance, 0.0, 1e-10);
}

NOVA_TEST(energy_balance_closes_mixed_source_and_neumann_flux) {
    const auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    ScalarDiffusion3DProblem p;
    p.mesh = box.mesh;
    p.diffusivity = [](double, double, double) { return 2.0; };
    p.source = [](double, double, double) { return 80.0; };
    p.cell_source.assign(p.mesh.elements.size(), 120.0);
    for (const auto& face : box.boundary_faces) {
        if (face_at_x(p.mesh, face, 0.0)) {
            p.dirichlet_faces.push_back({
                face, [](double, double, double) { return 400.0; }
            });
        } else if (face_at_x(p.mesh, face, 1.0)) {
            p.neumann_faces.push_back({
                face, [](double, double, double) { return 200.0; }
            });
        }
    }
    const auto solution = solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto budget = evaluate_steady_energy_balance_3d(p, solution.values);
    NOVA_REQUIRE_NEAR(budget.volumetric_generation, 200.0, 1e-9);
    NOVA_REQUIRE_NEAR(budget.neumann_outward, 200.0, 1e-9);
    NOVA_REQUIRE_NEAR(budget.dirichlet_outward, 0.0, 1e-8);
    NOVA_REQUIRE_NEAR(budget.total_outward, 200.0, 1e-8);
    NOVA_REQUIRE_NEAR(budget.imbalance, 0.0, 1e-8);
}

NOVA_TEST(energy_balance_recovers_convection_and_hot_wall_reaction) {
    const auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    ScalarDiffusion3DProblem p;
    p.mesh = box.mesh;
    p.diffusivity = [](double, double, double) { return 2.0; };
    p.source = [](double, double, double) { return 0.0; };
    for (const auto& face : box.boundary_faces) {
        if (face_at_x(p.mesh, face, 0.0)) {
            p.dirichlet_faces.push_back({
                face, [](double, double, double) { return 400.0; }
            });
        } else if (face_at_x(p.mesh, face, 1.0)) {
            p.robin_faces.push_back({
                face,
                [](double, double, double) { return 4.0; },
                [](double, double, double) { return 300.0; }
            });
        }
    }
    const auto solution = solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto budget = evaluate_steady_energy_balance_3d(p, solution.values);
    // 100 K across 1/k + 1/h = 0.75 K/W, unit face area.
    const double expected_rate = 100.0 / 0.75;
    NOVA_REQUIRE_NEAR(budget.robin_outward, expected_rate, 1e-8);
    NOVA_REQUIRE_NEAR(budget.dirichlet_outward, -expected_rate, 1e-8);
    NOVA_REQUIRE_NEAR(budget.volumetric_generation, 0.0, 1e-12);
    NOVA_REQUIRE_NEAR(budget.imbalance, 0.0, 1e-8);
}

NOVA_TEST(energy_balance_does_not_count_pure_neumann_gauge_as_physical_reaction) {
    auto p = unit_tetrahedron();
    p.source = [](double, double, double) { return 6.0; };
    p.neumann_faces.push_back({
        {0,2,3}, [](double, double, double) { return 2.0; }
    });
    const auto solution = solve_scalar_diffusion_3d(p);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto budget = evaluate_steady_energy_balance_3d(p, solution.values);
    NOVA_REQUIRE_NEAR(budget.volumetric_generation, 1.0, 1e-12);
    NOVA_REQUIRE_NEAR(budget.neumann_outward, 1.0, 1e-12);
    NOVA_REQUIRE_NEAR(budget.dirichlet_outward, 0.0, 1e-14);
    NOVA_REQUIRE_NEAR(budget.imbalance, 0.0, 1e-12);
}

NOVA_TEST(energy_balance_detects_nonconservative_supplied_free_node_field) {
    auto p = unit_tetrahedron();
    p.dirichlet_constraints = {{0, 0.0}};
    const std::vector<double> temperature{0.0, 1.0, 0.0, 0.0};
    const auto budget = evaluate_steady_energy_balance_3d(p, temperature);
    NOVA_REQUIRE_NEAR(budget.dirichlet_outward, 1.0/3.0, 1e-13);
    NOVA_REQUIRE_NEAR(budget.imbalance, 1.0/3.0, 1e-13);
}

NOVA_TEST(energy_balance_rejects_invalid_field_and_boundary_data) {
    auto p = unit_tetrahedron();
    const std::vector<double> values(4, 0.0);
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_steady_energy_balance_3d(p, std::vector<double>{0,0,0});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_steady_energy_balance_3d(
            p, std::vector<double>{0,0,0,std::numeric_limits<double>::quiet_NaN()});
    }));
    p.cell_source = {1.0, 2.0};
    NOVA_REQUIRE(rejects([&] { (void)evaluate_steady_energy_balance_3d(p, values); }));
    p.cell_source.clear();
    p.source = {};
    NOVA_REQUIRE(rejects([&] { (void)evaluate_steady_energy_balance_3d(p, values); }));
    p.source = [](double,double,double) { return 0.0; };
    p.robin_faces.push_back({
        {0,2,3},
        [](double,double,double) { return -1.0; },
        [](double,double,double) { return 0.0; }
    });
    NOVA_REQUIRE(rejects([&] { (void)evaluate_steady_energy_balance_3d(p, values); }));
}
