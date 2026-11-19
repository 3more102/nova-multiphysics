#include "test_framework.hpp"

#include "nova/fem/energy_balance_3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace nova::fem;

bool on_x(const TetrahedralMesh3D& mesh, const BoundaryFace3D& face, double x) {
    return std::all_of(face.begin(), face.end(),
                       [&](std::size_t node) { return mesh.nodes[node].x == x; });
}

ScalarDiffusion3DProblem wall(bool generated) {
    auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2, 2, 2);
    ScalarDiffusion3DProblem problem;
    problem.mesh = std::move(box.mesh);
    problem.source = [](double, double, double) { return 0.0; };
    for (const auto& cell : problem.mesh.elements) {
        double x = 0.0;
        for (auto node : cell) x += problem.mesh.nodes[node].x / 4;
        const bool left = x < 0.5;
        problem.cell_diffusivity.push_back(left ? 2.0 : 6.0);
        if (generated) problem.cell_source.push_back(left ? 120.0 : 360.0);
    }
    for (auto face : box.boundary_faces) {
        if (on_x(problem.mesh, face, 0)) {
            problem.dirichlet_faces.push_back({
                face, [](double, double, double) { return 400.0; }
            });
        } else if (on_x(problem.mesh, face, 1)) {
            problem.robin_faces.push_back({
                face,
                [](double, double, double) { return 12.0; },
                [](double, double, double) { return 300.0; }
            });
        }
    }
    return problem;
}

template<class F>
bool rejects(F&& f) {
    try { f(); }
    catch (const std::invalid_argument&) { return true; }
    return false;
}

} // namespace

NOVA_TEST(energy_balance_3d_exact_layered_wall_has_opposite_boundary_reactions) {
    const auto problem = wall(false);
    const auto solution = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto b = evaluate_diffusion_energy_balance_3d(problem, solution.values);
    NOVA_REQUIRE_NEAR(b.generated_heat, 0, 1e-11);
    NOVA_REQUIRE_NEAR(b.neumann_outward_heat, 0, 1e-11);
    NOVA_REQUIRE_NEAR(b.robin_outward_heat, 240, 2e-8);
    NOVA_REQUIRE_NEAR(b.constraint_outward_heat, -240, 2e-8);
    NOVA_REQUIRE_NEAR(b.signed_balance_defect, 0, 2e-8);
    NOVA_REQUIRE(b.max_free_equation_residual < 2e-8);
    NOVA_REQUIRE(b.constrained_nodes > 0);
    NOVA_REQUIRE(b.free_nodes > 0);
    NOVA_REQUIRE(b.constrained_nodes + b.free_nodes == problem.mesh.nodes.size());
}

NOVA_TEST(energy_balance_3d_regional_volumetric_heating_closes_global_budget) {
    const auto problem = wall(true);
    const auto solution = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto b = evaluate_diffusion_energy_balance_3d(problem, solution.values);
    NOVA_REQUIRE_NEAR(b.generated_heat, 240, 2e-9); // 120*.5+360*.5
    NOVA_REQUIRE_NEAR(b.neumann_outward_heat, 0, 2e-9);
    NOVA_REQUIRE_NEAR(b.robin_outward_heat, 384, 2e-7);
    NOVA_REQUIRE_NEAR(b.constraint_outward_heat, -144, 2e-7);
    NOVA_REQUIRE_NEAR(b.signed_balance_defect, 0, 2e-8);
    NOVA_REQUIRE(b.max_free_equation_residual < 2e-8);
}

NOVA_TEST(energy_balance_3d_neumann_and_spatial_source_follow_weak_form) {
    auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2, 2, 2);
    ScalarDiffusion3DProblem problem;
    problem.mesh = std::move(box.mesh);
    problem.diffusivity = [](double, double, double) { return 3.0; };
    problem.source = [](double x, double, double) { return 3.0 + 2.0*x; };
    for (const auto& face : box.boundary_faces) {
        if (on_x(problem.mesh, face, 0)) {
            problem.dirichlet_faces.push_back({
                face, [](double, double, double) { return 0.0; }
            });
        } else if (on_x(problem.mesh, face, 1)) {
            problem.neumann_faces.push_back({
                face, [](double, double, double) { return 4.0; }
            });
        }
    }
    const auto solution = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto b = evaluate_diffusion_energy_balance_3d(problem, solution.values);
    NOVA_REQUIRE_NEAR(b.generated_heat, 4.0, 1e-12);
    NOVA_REQUIRE_NEAR(b.neumann_outward_heat, 4.0, 1e-12);
    NOVA_REQUIRE_NEAR(b.robin_outward_heat, 0, 1e-12);
    NOVA_REQUIRE_NEAR(b.constraint_outward_heat, 0, 1e-8);
    NOVA_REQUIRE_NEAR(b.signed_balance_defect, 0, 1e-8);
    NOVA_REQUIRE(b.max_free_equation_residual < 1e-8);
}

NOVA_TEST(energy_balance_3d_detects_unconverged_field_even_if_global_heat_closes) {
    const auto problem = wall(false);
    auto solution = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(solution.linear_solver.converged);
    // Central free vertex. Keep all prescribed boundary values unchanged.
    auto i = std::find_if(problem.mesh.nodes.begin(), problem.mesh.nodes.end(),
                          [](const Point3D& p) {
                              return p.x == 0.5 && p.y == 0.5 && p.z == 0.5;
                          });
    NOVA_REQUIRE(i != problem.mesh.nodes.end());
    solution.values[static_cast<std::size_t>(i - problem.mesh.nodes.begin())] += 3.0;
    const auto b = evaluate_diffusion_energy_balance_3d(problem, solution.values);
    NOVA_REQUIRE(b.max_free_equation_residual > 0.01);
}


NOVA_TEST(energy_balance_3d_pure_neumann_zero_mean_gauge_is_not_a_reaction) {
    auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2, 2, 2);
    ScalarDiffusion3DProblem problem;
    problem.mesh = std::move(box.mesh);
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    const auto solution = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto b = evaluate_diffusion_energy_balance_3d(problem, solution.values);
    NOVA_REQUIRE(b.constrained_nodes == 0);
    NOVA_REQUIRE(b.free_nodes == problem.mesh.nodes.size());
    NOVA_REQUIRE_NEAR(b.constraint_outward_heat, 0.0, 1e-12);
    NOVA_REQUIRE_NEAR(b.signed_balance_defect, 0.0, 1e-12);
    NOVA_REQUIRE_NEAR(b.max_free_equation_residual, 0.0, 1e-12);
}

NOVA_TEST(energy_balance_3d_explicit_internal_pin_is_counted_as_a_reaction) {
    auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2, 2, 2);
    ScalarDiffusion3DProblem problem;
    problem.mesh = std::move(box.mesh);
    problem.diffusivity = [](double, double, double) { return 2.0; };
    problem.source = [](double, double, double) { return 10.0; };
    auto it = std::find_if(problem.mesh.nodes.begin(), problem.mesh.nodes.end(),
                           [](const Point3D& p) {
                               return p.x == 0.5 && p.y == 0.5 && p.z == 0.5;
                           });
    NOVA_REQUIRE(it != problem.mesh.nodes.end());
    const auto center = static_cast<std::size_t>(it - problem.mesh.nodes.begin());
    problem.dirichlet_constraints.push_back({center, 0.0});
    const auto solution = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto b = evaluate_diffusion_energy_balance_3d(problem, solution.values);
    NOVA_REQUIRE(b.constrained_nodes == 1);
    NOVA_REQUIRE_NEAR(b.generated_heat, 10.0, 1e-10);
    NOVA_REQUIRE_NEAR(b.constraint_outward_heat, 10.0, 1e-8);
    NOVA_REQUIRE_NEAR(b.signed_balance_defect, 0.0, 1e-8);
    NOVA_REQUIRE(b.max_free_equation_residual < 1e-8);
}


NOVA_TEST(energy_balance_3d_anisotropic_spd_materials_preserve_heating_budget) {
    auto problem = wall(true);
    for (double longitudinal : problem.cell_diffusivity) {
        problem.cell_diffusivity_tensors.push_back({
            longitudinal, 4.0, 5.0, 0.0, 0.0, 0.0
        });
    }
    problem.cell_diffusivity.clear();
    const auto solution = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(solution.linear_solver.converged);
    const auto b = evaluate_diffusion_energy_balance_3d(problem, solution.values);
    NOVA_REQUIRE_NEAR(b.generated_heat, 240.0, 2e-9);
    NOVA_REQUIRE_NEAR(b.robin_outward_heat, 384.0, 2e-7);
    NOVA_REQUIRE_NEAR(b.constraint_outward_heat, -144.0, 2e-7);
    NOVA_REQUIRE_NEAR(b.signed_balance_defect, 0.0, 2e-8);
    NOVA_REQUIRE(b.max_free_equation_residual < 2e-8);
}

NOVA_TEST(energy_balance_3d_validates_sources_mesh_boundary_and_constraints) {
    auto problem = wall(false);
    auto solution = solve_scalar_diffusion_3d(problem);
    const auto field = solution.values;
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate_diffusion_energy_balance_3d(
            problem, std::vector<double>{1.0, 2.0});
    }));
    auto modified = field;
    modified[0] = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] { (void)evaluate_diffusion_energy_balance_3d(problem, modified); }));
    problem.source = {};
    NOVA_REQUIRE(rejects([&] { (void)evaluate_diffusion_energy_balance_3d(problem, field); }));
    problem.source = [](double, double, double) { return 0.0; };
    problem.cell_source = {1.0};
    NOVA_REQUIRE(rejects([&] { (void)evaluate_diffusion_energy_balance_3d(problem, field); }));
    problem.cell_source.clear();
    problem.source = [](double, double, double) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    NOVA_REQUIRE(rejects([&] { (void)evaluate_diffusion_energy_balance_3d(problem, field); }));
    problem.source = [](double, double, double) { return 0.0; };
    problem.robin_faces.push_back(problem.robin_faces.front());
    NOVA_REQUIRE(rejects([&] { (void)evaluate_diffusion_energy_balance_3d(problem, field); }));
    problem.robin_faces.pop_back();
    auto fixed = field;
    auto first = problem.dirichlet_faces.front().face.front();
    fixed[first] += 1.0;
    NOVA_REQUIRE(rejects([&] { (void)evaluate_diffusion_energy_balance_3d(problem, fixed); }));
    problem.diffusivity = [](double, double, double) { return -1.0; };
    problem.cell_diffusivity.clear();
    NOVA_REQUIRE(rejects([&] { (void)evaluate_diffusion_energy_balance_3d(problem, field); }));
}
