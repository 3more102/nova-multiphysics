#include "test_framework.hpp"

#include "nova/fem/adaptive_diffusion_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace {

[[nodiscard]] bool is_x_plate(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::BoundaryFace3D& face, double x
) {
    return mesh.nodes[face[0]].x == x &&
           mesh.nodes[face[1]].x == x &&
           mesh.nodes[face[2]].x == x;
}

template <class F>
[[nodiscard]] bool rejects(F&& function) {
    try {
        function();
    } catch (const std::exception&) {
        return true;
    }
    return false;
}

[[nodiscard]] nova::fem::ScalarDiffusion3DProblem make_wall(
    const nova::fem::StructuredBoxMesh3D& box,
    double source
) {
    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    problem.source = [source](double, double, double) { return source; };
    for (const auto& face : box.boundary_faces) {
        if (is_x_plate(box.mesh, face, 0.0) ||
            is_x_plate(box.mesh, face, 1.0)) {
            problem.dirichlet_faces.push_back({
                face, [](double, double, double) { return 0.0; }});
        }
    }
    return problem;
}

// Independently integrate k * |grad(u_exact)-grad(u_h)|^2 on each
// tetrahedron for u_exact=x(1-x), k=2, -div(k grad u)=4. Because
// grad(u_h) is cellwise constant and grad(u_exact) is affine in x,
// first and second tetrahedral x moments give the exact integral.
[[nodiscard]] double manufactured_wall_energy_error_squared(
    const nova::fem::ScalarDiffusion3DProblem& problem,
    const std::vector<double>& values
) {
    const auto fluxes = nova::fem::evaluate_p1_tetrahedral_flux_3d(
        problem, values);
    const auto quality = nova::fem::evaluate_mesh_quality_3d(problem.mesh);
    double error_squared = 0.0;
    for (std::size_t i = 0; i < problem.mesh.elements.size(); ++i) {
        double sum_x = 0.0;
        double sum_x_squared = 0.0;
        for (const auto node : problem.mesh.elements[i]) {
            const double x = problem.mesh.nodes[node].x;
            sum_x += x;
            sum_x_squared += x * x;
        }
        const double mean_x = sum_x / 4.0;
        const double mean_x_squared =
            (sum_x * sum_x + sum_x_squared) / 20.0;
        const auto& gradient = fluxes[i].gradient;
        const double a = 1.0 - gradient.x;
        const double mean_squared_gradient_error =
            a * a - 4.0 * a * mean_x + 4.0 * mean_x_squared +
            gradient.y * gradient.y + gradient.z * gradient.z;
        error_squared += problem.cell_diffusivity[i] *
            quality.elements[i].volume * mean_squared_gradient_error;
    }
    return error_squared;
}

} // namespace

NOVA_TEST(adaptive_diffusion_3d_solves_marks_refines_transfers_and_resolves) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const auto problem = make_wall(box, 4.0);
    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 0.5);

    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.adaptation.refinement.mesh.nodes.size() >
                 problem.mesh.nodes.size());
    NOVA_REQUIRE(result.adaptation.refinement.mesh.elements.size() >
                 problem.mesh.elements.size());
    result.adaptation.refinement.mesh.validate();
    NOVA_REQUIRE(result.refined_problem.mesh.elements.size() ==
                 result.adaptation.refinement.mesh.elements.size());
    NOVA_REQUIRE(result.refined_problem.cell_diffusivity.size() ==
                 result.refined_problem.mesh.elements.size());
    NOVA_REQUIRE(result.coarse_indicators.total_error > 0.0);
    NOVA_REQUIRE(std::isfinite(result.refined_indicators.total_error));
    NOVA_REQUIRE(result.refined_indicators.element_errors.size() ==
                 result.refined_problem.mesh.elements.size());
    for (const auto k : result.refined_problem.cell_diffusivity) {
        NOVA_REQUIRE_NEAR(k, 2.0, 1.0e-14);
    }
    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.transferred_p1_values[i],
            result.coarse_solution.values[i], 1.0e-14);
    }
    for (std::size_t i = 0; i < result.adaptation.refinement.midpoint_parent_edges.size(); ++i) {
        const auto& edge = result.adaptation.refinement.midpoint_parent_edges[i];
        NOVA_REQUIRE_NEAR(
            result.transferred_p1_values[problem.mesh.nodes.size() + i],
            0.5 * (result.coarse_solution.values[edge[0]] +
                   result.coarse_solution.values[edge[1]]), 1.0e-14);
    }
    NOVA_REQUIRE(result.refined_problem.dirichlet_faces.size() >
                 problem.dirichlet_faces.size());
    // For 2*u'' = -4, exact u=x(1-x); high-level sanity and essential BCs.
    for (std::size_t i = 0; i < result.refined_problem.mesh.nodes.size(); ++i) {
        const auto& p = result.refined_problem.mesh.nodes[i];
        if (p.x == 0.0 || p.x == 1.0) {
            NOVA_REQUIRE_NEAR(result.refined_solution.values[i], 0.0, 1.0e-8);
        }
    }
}

NOVA_TEST(adaptive_diffusion_3d_reduces_exact_manufactured_energy_error) {
    // The source-loaded wall has an exact continuum solution u=x(1-x).
    // Independently measure its energy-seminorm error, rather than relying
    // solely on the residual indicator's own self-reported decrease.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const auto problem = make_wall(box, 4.0);
    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 0.65);

    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    const double coarse_error_squared =
        manufactured_wall_energy_error_squared(
            problem, result.coarse_solution.values);
    const double fine_error_squared =
        manufactured_wall_energy_error_squared(
            result.refined_problem, result.refined_solution.values);
    NOVA_REQUIRE(std::isfinite(coarse_error_squared));
    NOVA_REQUIRE(std::isfinite(fine_error_squared));
    NOVA_REQUIRE(coarse_error_squared > 1.0e-6);
    NOVA_REQUIRE(fine_error_squared >= 0.0);
    NOVA_REQUIRE(fine_error_squared < coarse_error_squared - 1.0e-9);
    NOVA_REQUIRE_NEAR(
        nova::fem::evaluate_mesh_quality_3d(
            result.refined_problem.mesh).total_volume,
        1.0, 1.0e-12);
}

NOVA_TEST(adaptive_diffusion_3d_zero_indicator_returns_identity_without_resolve) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    auto problem = make_wall(box, 0.0);
    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 0.5);
    NOVA_REQUIRE(result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.adaptation.refinement.mesh.nodes.size() ==
                 problem.mesh.nodes.size());
    NOVA_REQUIRE(result.adaptation.refinement.mesh.elements.size() ==
                 problem.mesh.elements.size());
    NOVA_REQUIRE(result.coarse_indicators.total_error == 0.0);
    NOVA_REQUIRE(result.refined_indicators.total_error == 0.0);
    NOVA_REQUIRE(result.refined_solution.values ==
                 result.coarse_solution.values);
    NOVA_REQUIRE(result.transferred_p1_values ==
                 result.coarse_solution.values);
}

NOVA_TEST(adaptive_diffusion_3d_transfers_prescribed_neumann_faces) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    auto problem = make_wall(box, 4.0);
    for (const auto& face : box.boundary_faces) {
        const auto& mesh = box.mesh;
        if (mesh.nodes[face[0]].y == 0.0 &&
            mesh.nodes[face[1]].y == 0.0 &&
            mesh.nodes[face[2]].y == 0.0) {
            problem.neumann_faces.push_back({
                face, [](double, double, double) { return 1.0; }});
        }
    }
    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 0.5);
    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.refined_problem.neumann_faces.size() >
                 problem.neumann_faces.size());
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    NOVA_REQUIRE(std::isfinite(result.refined_indicators.total_error));
    for (const auto& face : result.refined_problem.neumann_faces) {
        NOVA_REQUIRE_NEAR(
            face.outward_flux(0.0, 0.0, 0.0), 1.0, 1.0e-14);
    }
}

NOVA_TEST(adaptive_diffusion_3d_rejects_unrepresented_physics) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    auto problem = make_wall(box, 1.0);
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::solve_adaptive_diffusion_step_3d(problem, 0.0);
    }));
    auto modified = problem;
    modified.cell_diffusivity.clear();
    modified.diffusivity = [](double, double, double) { return 2.0; };
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::solve_adaptive_diffusion_step_3d(modified, 0.5);
    }));
    modified = problem;
    modified.cell_source.assign(box.mesh.elements.size(), 1.0);
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::solve_adaptive_diffusion_step_3d(modified, 0.5);
    }));
    modified = problem;
    modified.robin_faces.push_back({
        box.boundary_faces.front(),
        [](double, double, double) { return 1.0; },
        [](double, double, double) { return 0.0; }});
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::solve_adaptive_diffusion_step_3d(modified, 0.5);
    }));
}


NOVA_TEST(adaptive_diffusion_3d_re_evaluates_spatial_boundary_values_after_refinement) {
    // Manufactured solution u=x(1-x)+y^2/4, isotropic K=2:
    // -div(K grad u)=3; on x=0,1 u=y^2/4; on y=1 outward
    // flux=-1; y=0 and z-normal surfaces are insulated.
    // A fine boundary midpoint must use the callback's actual value,
    // rather than interpolation of the coarse boundary nodal values.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    problem.source = [](double, double, double) { return 3.0; };
    for (const auto& face : box.boundary_faces) {
        if (is_x_plate(box.mesh, face, 0.0) ||
            is_x_plate(box.mesh, face, 1.0)) {
            problem.dirichlet_faces.push_back({
                face, [](double, double y, double) { return 0.25 * y * y; }});
        } else if (box.mesh.nodes[face[0]].y == 1.0 &&
                   box.mesh.nodes[face[1]].y == 1.0 &&
                   box.mesh.nodes[face[2]].y == 1.0) {
            problem.neumann_faces.push_back({
                face, [](double, double, double) { return -1.0; }});
        }
    }

    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 1.0);
    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.refined_problem.dirichlet_faces.size() >
                 problem.dirichlet_faces.size());
    NOVA_REQUIRE(result.refined_problem.neumann_faces.size() >
                 problem.neumann_faces.size());

    bool new_dirichlet_vertex = false;
    bool coarse_interpolation_differs = false;
    for (std::size_t i = problem.mesh.nodes.size();
         i < result.refined_problem.mesh.nodes.size(); ++i) {
        const auto& p = result.refined_problem.mesh.nodes[i];
        if (p.x != 0.0 && p.x != 1.0) continue;
        new_dirichlet_vertex = true;
        const double exact_boundary = 0.25 * p.y * p.y;
        NOVA_REQUIRE_NEAR(result.refined_solution.values[i],
                          exact_boundary, 1.0e-8);
        if (std::abs(result.transferred_p1_values[i] -
                     exact_boundary) > 1.0e-5) {
            coarse_interpolation_differs = true;
        }
    }
    NOVA_REQUIRE(new_dirichlet_vertex);
    NOVA_REQUIRE(coarse_interpolation_differs);
    for (const auto& condition : result.refined_problem.neumann_faces) {
        NOVA_REQUIRE_NEAR(condition.outward_flux(0.5, 1.0, 0.5),
                          -1.0, 1.0e-14);
    }
}
