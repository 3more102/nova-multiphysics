#include "test_framework.hpp"

#include "nova/fem/adaptive_diffusion_3d.hpp"
#include "nova/fem/p1_flux_3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
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

NOVA_TEST(adaptive_diffusion_3d_rejects_invalid_parameter_and_material_selection) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    auto problem = make_wall(box, 1.0);
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::solve_adaptive_diffusion_step_3d(problem, 0.0);
    }));
    auto modified = problem;
    // Keep the cellwise coefficient while adding a spatial coefficient: the
    // base solver must reject selecting two material representations.
    modified.diffusivity = [](double, double, double) { return 2.0; };
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


NOVA_TEST(adaptive_diffusion_3d_heterogeneous_cell_materials_and_natural_face_coverage) {
    // Exercise both supported heterogeneous cellwise k and the implicit
    // zero-Neumann faces that would otherwise be easy to omit from eta_K.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    auto problem = make_wall(box, 4.0);
    for (std::size_t i = 0; i < problem.cell_diffusivity.size(); ++i) {
        problem.cell_diffusivity[i] = 1.0 + static_cast<double>(i % 4U);
    }

    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 0.6);
    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());

    const auto& refinement = result.adaptation.refinement;
    NOVA_REQUIRE(result.refined_problem.cell_diffusivity.size() ==
                 refinement.parent_elements.size());
    for (std::size_t child = 0; child < refinement.parent_elements.size(); ++child) {
        const std::size_t parent = refinement.parent_elements[child];
        NOVA_REQUIRE(parent < problem.cell_diffusivity.size());
        NOVA_REQUIRE(result.refined_problem.cell_diffusivity[child] ==
                     problem.cell_diffusivity[parent]);
    }

    // Each exterior face is either explicitly Dirichlet or is reported
    // as natural, including omitted insulated faces with zero flux.
    const auto verify_natural_coverage = [](
        const nova::fem::ScalarDiffusion3DProblem& model,
        const nova::fem::DiffusionResidualIndicators3D& indicators
    ) {
        const auto exterior = nova::fem::extract_boundary_faces_3d(model.mesh);
        NOVA_REQUIRE(indicators.boundary_residual_errors.size() ==
                     model.mesh.elements.size());
        std::set<nova::fem::BoundaryFace3D> exterior_keys;
        for (auto face : exterior) {
            std::sort(face.begin(), face.end());
            NOVA_REQUIRE(exterior_keys.insert(face).second);
        }
        std::set<nova::fem::BoundaryFace3D> classified;
        for (const auto& condition : model.dirichlet_faces) {
            auto face = condition.face;
            std::sort(face.begin(), face.end());
            NOVA_REQUIRE(exterior_keys.contains(face));
            NOVA_REQUIRE(classified.insert(face).second);
        }
        for (const auto& report : indicators.boundary_residual_faces) {
            NOVA_REQUIRE(exterior_keys.contains(report.nodes));
            NOVA_REQUIRE(classified.insert(report.nodes).second);
            const auto& a = model.mesh.nodes[report.nodes[0]];
            const auto& b = model.mesh.nodes[report.nodes[1]];
            const auto& c = model.mesh.nodes[report.nodes[2]];
            // The x-end plates are Dirichlet; no plate belongs to the
            // natural-boundary residual estimator.
            const bool x_plate =
                (a.x == 0.0 && b.x == 0.0 && c.x == 0.0) ||
                (a.x == 1.0 && b.x == 1.0 && c.x == 1.0);
            NOVA_REQUIRE(!x_plate);
            NOVA_REQUIRE(std::isfinite(report.mismatch_l2));
            NOVA_REQUIRE(report.mismatch_l2 >= 0.0);
        }
        NOVA_REQUIRE(classified.size() == exterior_keys.size());
    };
    verify_natural_coverage(problem, result.coarse_indicators);
    verify_natural_coverage(result.refined_problem, result.refined_indicators);
}


NOVA_TEST(adaptive_diffusion_3d_vertex_only_affine_dirichlet_is_not_neumann) {
    // With every boundary vertex fixed to u=x, each boundary triangle has
    // an essential P1 trace. It must not contribute an insulated-face
    // flux residual (which would cause spurious adaptive refinement).
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    for (std::size_t i = 0; i < box.mesh.nodes.size(); ++i) {
        problem.dirichlet_constraints.push_back({
            i, box.mesh.nodes[i].x});
    }
    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 0.5);
    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.coarse_indicators.boundary_residual_faces.empty());
    NOVA_REQUIRE_NEAR(result.coarse_indicators.total_error, 0.0, 1.0e-10);
    NOVA_REQUIRE(result.adaptation.marked_elements.empty());
}

NOVA_TEST(adaptive_diffusion_3d_preserves_vertex_only_dirichlet_faces) {
    // A positive source triggers refinement. All boundary vertices are
    // prescribed u=y, without face callbacks; each new midpoint on those
    // fully fixed boundary faces must retain the P1 interpolated trace.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    problem.source = [](double, double, double) { return 4.0; };
    for (std::size_t i = 0; i < box.mesh.nodes.size(); ++i) {
        problem.dirichlet_constraints.push_back({
            i, box.mesh.nodes[i].y});
    }
    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(
        problem, 1.0);
    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.refined_problem.dirichlet_constraints.size() >
                 problem.dirichlet_constraints.size());
    NOVA_REQUIRE(result.refined_indicators.boundary_residual_faces.empty());
    for (const auto& face : result.adaptation.refinement.boundary_faces) {
        for (const auto node : face) {
            const auto& p = result.refined_problem.mesh.nodes[node];
            NOVA_REQUIRE_NEAR(result.refined_solution.values[node],
                              p.y, 1.0e-8);
        }
    }
}


NOVA_TEST(adaptive_diffusion_3d_spatial_tensor_reconstructs_exact_off_diagonal_flux) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const nova::fem::SymmetricTensor3D tensor{4.0, 3.0, 2.0, 0.5, 0.1, 0.2};
    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity_tensor = [tensor](double, double, double) { return tensor; };
    for (std::size_t node = 0; node < box.mesh.nodes.size(); ++node) {
        const auto& point = box.mesh.nodes[node];
        problem.dirichlet_constraints.push_back({
            node, point.x + 2.0 * point.y - 3.0 * point.z});
    }

    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(problem, 0.5);
    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    NOVA_REQUIRE(result.coarse_indicators.total_error < 1.0e-10);

    const auto fields = nova::fem::evaluate_p1_tetrahedral_flux_3d(
        problem, result.coarse_solution.values);
    NOVA_REQUIRE(fields.size() == box.mesh.elements.size());
    for (const auto& field : fields) {
        NOVA_REQUIRE_NEAR(field.gradient.x, 1.0, 1.0e-11);
        NOVA_REQUIRE_NEAR(field.gradient.y, 2.0, 1.0e-11);
        NOVA_REQUIRE_NEAR(field.gradient.z, -3.0, 1.0e-11);
        NOVA_REQUIRE_NEAR(field.flux.x, -4.7, 1.0e-10);
        NOVA_REQUIRE_NEAR(field.flux.y, -5.9, 1.0e-10);
        NOVA_REQUIRE_NEAR(field.flux.z, 5.5, 1.0e-10);
    }
}

NOVA_TEST(adaptive_diffusion_3d_transfers_cell_tensors_and_robin_neumann_faces) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const nova::fem::SymmetricTensor3D tensor{4.0, 3.0, 2.0, 0.5, 0.1, 0.2};
    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_diffusivity_tensors.assign(box.mesh.elements.size(), tensor);
    problem.source = [](double, double, double) { return 4.0; };

    const auto on_y = [&](const nova::fem::BoundaryFace3D& face, double y) {
        return box.mesh.nodes[face[0]].y == y &&
               box.mesh.nodes[face[1]].y == y &&
               box.mesh.nodes[face[2]].y == y;
    };
    for (const auto& face : box.boundary_faces) {
        if (is_x_plate(box.mesh, face, 0.0) || is_x_plate(box.mesh, face, 1.0)) {
            problem.dirichlet_faces.push_back({
                face, [](double, double, double) { return 0.0; }});
        } else if (on_y(face, 1.0)) {
            problem.robin_faces.push_back({
                face,
                [](double, double, double) { return 1.5; },
                [](double, double, double) { return 0.25; }});
        } else if (on_y(face, 0.0)) {
            problem.neumann_faces.push_back({
                face, [](double, double, double) { return 0.5; }});
        }
    }

    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(problem, 1.0);
    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.refined_problem.cell_diffusivity_tensors.size() ==
                 result.refined_problem.mesh.elements.size());
    NOVA_REQUIRE(result.refined_problem.robin_faces.size() > problem.robin_faces.size());
    NOVA_REQUIRE(result.refined_problem.neumann_faces.size() > problem.neumann_faces.size());

    const auto& refinement = result.adaptation.refinement;
    NOVA_REQUIRE(refinement.parent_elements.size() ==
                 result.refined_problem.cell_diffusivity_tensors.size());
    for (std::size_t child = 0; child < refinement.parent_elements.size(); ++child) {
        const auto parent = refinement.parent_elements[child];
        NOVA_REQUIRE(parent < problem.cell_diffusivity_tensors.size());
        const auto& actual = result.refined_problem.cell_diffusivity_tensors[child];
        const auto& expected = problem.cell_diffusivity_tensors[parent];
        NOVA_REQUIRE_NEAR(actual.xx, expected.xx, 1.0e-14);
        NOVA_REQUIRE_NEAR(actual.yy, expected.yy, 1.0e-14);
        NOVA_REQUIRE_NEAR(actual.zz, expected.zz, 1.0e-14);
        NOVA_REQUIRE_NEAR(actual.xy, expected.xy, 1.0e-14);
        NOVA_REQUIRE_NEAR(actual.xz, expected.xz, 1.0e-14);
        NOVA_REQUIRE_NEAR(actual.yz, expected.yz, 1.0e-14);
    }

    bool saw_robin_residual = false;
    bool saw_neumann_residual = false;
    for (const auto& face : result.refined_problem.robin_faces) {
        NOVA_REQUIRE_NEAR(face.transfer_coefficient(0.25, 1.0, 0.5), 1.5, 1.0e-14);
        NOVA_REQUIRE_NEAR(face.ambient_value(0.25, 1.0, 0.5), 0.25, 1.0e-14);
    }
    for (const auto& face : result.refined_indicators.boundary_residual_faces) {
        saw_robin_residual = saw_robin_residual ||
            face.kind == nova::fem::BoundaryResidualKind3D::Robin;
        saw_neumann_residual = saw_neumann_residual ||
            face.kind == nova::fem::BoundaryResidualKind3D::Neumann;
    }
    NOVA_REQUIRE(saw_robin_residual);
    NOVA_REQUIRE(saw_neumann_residual);
    NOVA_REQUIRE(std::isfinite(result.refined_indicators.total_error));
}


NOVA_TEST(adaptive_diffusion_3d_transfers_additive_cell_source_to_children) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    nova::fem::ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.cell_diffusivity.assign(box.mesh.elements.size(), 2.0);
    problem.source = [](double, double, double) { return 1.0; };
    problem.cell_source.resize(box.mesh.elements.size());
    for (std::size_t cell = 0; cell < problem.cell_source.size(); ++cell) {
        problem.cell_source[cell] = 0.25 * static_cast<double>(cell + 1U);
    }
    for (const auto& face : box.boundary_faces) {
        if (is_x_plate(box.mesh, face, 0.0) || is_x_plate(box.mesh, face, 1.0)) {
            problem.dirichlet_faces.push_back({
                face, [](double, double, double) { return 0.0; }});
        }
    }

    const auto result = nova::fem::solve_adaptive_diffusion_step_3d(problem, 1.0);
    NOVA_REQUIRE(result.coarse_solution.linear_solver.converged);
    NOVA_REQUIRE(result.refined_solution.linear_solver.converged);
    NOVA_REQUIRE(!result.adaptation.marked_elements.empty());
    NOVA_REQUIRE(result.refined_problem.cell_source.size() ==
                 result.refined_problem.mesh.elements.size());
    NOVA_REQUIRE_NEAR(result.refined_problem.source(0.2, 0.3, 0.4), 1.0, 1.0e-14);

    // Independently verify h_K * ||f_spatial + f_cell||_L2(K) for the
    // piecewise-constant cell source plus the constant spatial source.
    const auto coarse_quality = nova::fem::evaluate_mesh_quality_3d(problem.mesh);
    NOVA_REQUIRE(result.coarse_indicators.cell_residual_errors.size() ==
                 problem.mesh.elements.size());
    for (std::size_t cell = 0; cell < problem.mesh.elements.size(); ++cell) {
        const auto& tetra = problem.mesh.elements[cell];
        double diameter = 0.0;
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) {
                const auto& a = problem.mesh.nodes[tetra[i]];
                const auto& b = problem.mesh.nodes[tetra[j]];
                diameter = std::max(diameter, std::hypot(
                    b.x - a.x, b.y - a.y, b.z - a.z));
            }
        }
        const double source_value = 1.0 + problem.cell_source[cell];
        const double expected = diameter *
            std::sqrt(coarse_quality.elements[cell].volume) * source_value;
        NOVA_REQUIRE_NEAR(
            result.coarse_indicators.cell_residual_errors[cell], expected, 1.0e-10);
    }

    const auto& refinement = result.adaptation.refinement;
    NOVA_REQUIRE(refinement.parent_elements.size() ==
                 result.refined_problem.cell_source.size());
    for (std::size_t child = 0; child < refinement.parent_elements.size(); ++child) {
        const auto parent = refinement.parent_elements[child];
        NOVA_REQUIRE(parent < problem.cell_source.size());
        NOVA_REQUIRE_NEAR(result.refined_problem.cell_source[child],
                          problem.cell_source[parent], 1.0e-14);
    }
    NOVA_REQUIRE(std::isfinite(result.coarse_indicators.total_error));
    NOVA_REQUIRE(std::isfinite(result.refined_indicators.total_error));
}


NOVA_TEST(adaptive_diffusion_3d_multicycle_driver_stops_when_indicator_tolerance_is_met) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const auto problem = make_wall(box, 0.0);
    nova::fem::AdaptiveDiffusionOptions3D options;
    options.max_cycles = 4U;
    options.absolute_tolerance = 1.0e-12;
    options.relative_tolerance = 1.0e-3;

    const auto result = nova::fem::solve_adaptive_diffusion_3d(problem, options);
    NOVA_REQUIRE(result.tolerance_satisfied);
    NOVA_REQUIRE(result.stop_reason ==
                 nova::fem::AdaptiveDiffusionStopReason3D::ToleranceSatisfied);
    NOVA_REQUIRE(result.cycles.size() == 1U);
    NOVA_REQUIRE(result.cycles.front().marked_elements == 0U);
    NOVA_REQUIRE(result.final_problem.mesh.elements.size() ==
                 problem.mesh.elements.size());
    NOVA_REQUIRE(result.final_solution.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.final_indicators.total_error, 0.0, 1.0e-12);

    options.max_cycles = 0U;
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::solve_adaptive_diffusion_3d(problem, options);
    }));
}

NOVA_TEST(adaptive_diffusion_3d_multicycle_driver_respects_cycle_cap_and_reports_progress) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const auto problem = make_wall(box, 4.0);
    nova::fem::AdaptiveDiffusionOptions3D options;
    options.dorfler_theta = 0.6;
    options.absolute_tolerance = 0.0;
    options.relative_tolerance = 0.0;
    options.max_cycles = 2U;

    const auto result = nova::fem::solve_adaptive_diffusion_3d(problem, options);
    NOVA_REQUIRE(!result.tolerance_satisfied);
    NOVA_REQUIRE(result.stop_reason ==
                 nova::fem::AdaptiveDiffusionStopReason3D::MaximumCyclesReached);
    NOVA_REQUIRE(result.cycles.size() == options.max_cycles);
    NOVA_REQUIRE(result.final_solution.linear_solver.converged);
    NOVA_REQUIRE(result.final_problem.mesh.elements.size() >
                 problem.mesh.elements.size());
    for (const auto& cycle : result.cycles) {
        NOVA_REQUIRE(cycle.marked_elements > 0U);
        NOVA_REQUIRE(cycle.refined_elements > cycle.coarse_elements);
        NOVA_REQUIRE(cycle.refined_nodes > cycle.coarse_nodes);
        NOVA_REQUIRE(cycle.coarse_indicator > 0.0);
        NOVA_REQUIRE(cycle.refined_indicator > 0.0);
    }
    NOVA_REQUIRE(result.final_problem.mesh.elements.size() ==
                 result.cycles.back().refined_elements);
}
