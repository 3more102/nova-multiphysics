#include "test_framework.hpp"

#include "nova/fem/adaptive_diffusion_3d.hpp"
#include "nova/fem/error_norms_3d.hpp"

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


NOVA_TEST(adaptive_diffusion_3d_manufactured_wall_reduces_independent_l2_and_h1_error) {
    // -div(2 grad u)=4 on the unit cube with u=0 at x=0,1 and
    // zero flux elsewhere has the exact solution u=x(1-x).
    // theta=1 marks all nonzero indicators for this smooth benchmark;
    // this verifies one global adaptive step, not monotonicity in general.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const auto problem = make_wall(box, 4.0);
    const nova::fem::ExactScalarField3D exact{
        [](double x, double, double) { return x * (1.0 - x); },
        [](double x, double, double) { return 1.0 - 2.0 * x; },
        [](double, double, double) { return 0.0; },
        [](double, double, double) { return 0.0; }
    };

    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 1.0);
    NOVA_REQUIRE(result.adaptation.marked_elements.size() ==
                 problem.mesh.elements.size());
    const auto coarse = nova::fem::compute_error_norms_3d(
        problem.mesh, result.coarse_solution.values, exact);
    const auto refined = nova::fem::compute_error_norms_3d(
        result.refined_problem.mesh, result.refined_solution.values, exact);
    NOVA_REQUIRE(coarse.l2 > 0.0);
    NOVA_REQUIRE(coarse.h1_seminorm > 0.0);
    NOVA_REQUIRE(refined.l2 < 0.9 * coarse.l2);
    NOVA_REQUIRE(refined.h1_seminorm < 0.9 * coarse.h1_seminorm);
}

NOVA_TEST(adaptive_diffusion_3d_reevaluates_nonconstant_dirichlet_faces_on_new_nodes) {
    // The exact u=x^2+y^2+z^2 solves -Delta u=-6. A curved
    // prescribed boundary field must be sampled at NEW P1 boundary
    // nodes, rather than interpolated from old boundary vertices.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const auto exact = [](double x, double y, double z) {
        return x*x + y*y + z*z;
    };
    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_diffusivity.assign(box.mesh.elements.size(), 1.0);
    problem.source = [](double, double, double) { return -6.0; };
    for (const auto& face : box.boundary_faces) {
        problem.dirichlet_faces.push_back({face, exact});
    }

    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 0.5);
    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.refined_problem.dirichlet_faces.size() >
                 problem.dirichlet_faces.size());
    bool saw_new_boundary_vertex = false;
    for (const auto& condition : result.refined_problem.dirichlet_faces) {
        for (const auto node : condition.face) {
            const auto& p = result.refined_problem.mesh.nodes[node];
            NOVA_REQUIRE_NEAR(condition.value(p.x, p.y, p.z),
                              exact(p.x, p.y, p.z), 1.0e-14);
            NOVA_REQUIRE_NEAR(result.refined_solution.values[node],
                              exact(p.x, p.y, p.z), 1.0e-8);
            if (node >= problem.mesh.nodes.size()) {
                saw_new_boundary_vertex = true;
            }
        }
    }
    NOVA_REQUIRE(saw_new_boundary_vertex);
}
