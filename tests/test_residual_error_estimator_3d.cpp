#include "test_framework.hpp"

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/error_norms_3d.hpp"
#include "nova/fem/mesh_refinement_3d.hpp"
#include "nova/fem/h1_space_3d.hpp"
#include "nova/fem/residual_error_estimator_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using nova::fem::H1FunctionSpace3D;
using nova::fem::TetrahedralMesh3D;

template <typename Function>
[[nodiscard]] bool rejects(Function&& function) {
    try {
        function();
    } catch (const std::exception&) {
        return true;
    }
    return false;
}

[[nodiscard]] TetrahedralMesh3D two_tetrahedra_sharing_face() {
    TetrahedralMesh3D mesh;
    mesh.nodes = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 1.0}
    };
    mesh.elements = {{0U, 1U, 2U, 3U}, {1U, 2U, 3U, 4U}};
    return mesh;
}

} // namespace

NOVA_TEST(diffusion_residual_estimator_3d_is_zero_for_global_affine_harmonic_field) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const H1FunctionSpace3D space(box.mesh);
    const auto values = space.interpolate(
        [](double x, double y, double z) {
            return 1.0 + 2.0 * x - 3.0 * y + 4.0 * z;
        });
    const std::vector<double> conductivity(box.mesh.elements.size(), 1.0);
    const auto result = nova::fem::estimate_diffusion_residual_error_3d(
        box.mesh, values, conductivity,
        [](double, double, double) { return 0.0; });
    NOVA_REQUIRE(result.element_errors.size() == box.mesh.elements.size());
    NOVA_REQUIRE(result.total_error < 1.0e-12);
    NOVA_REQUIRE(result.flux_jumps.faces.size() > 0U);
}

NOVA_TEST(diffusion_residual_estimator_3d_integrates_constant_cell_source) {
    TetrahedralMesh3D mesh;
    mesh.nodes = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };
    mesh.elements = {{0U, 1U, 2U, 3U}};
    const H1FunctionSpace3D space(mesh);
    const auto values = space.interpolate(
        [](double, double, double) { return 0.0; });
    const auto result = nova::fem::estimate_diffusion_residual_error_3d(
        mesh, values, std::vector<double>{1.0},
        [](double, double, double) { return 1.0; });
    const double exact = std::sqrt(1.0 / 3.0);
    NOVA_REQUIRE_NEAR(result.cell_residual_errors[0], exact, 1.0e-13);
    NOVA_REQUIRE_NEAR(result.element_errors[0], exact, 1.0e-13);
    NOVA_REQUIRE_NEAR(result.total_error, exact, 1.0e-13);
}

NOVA_TEST(diffusion_residual_estimator_3d_detects_interior_flux_jump) {
    const auto mesh = two_tetrahedra_sharing_face();
    const H1FunctionSpace3D space(mesh);
    std::vector<double> values(space.dimension(), 0.0);
    values[4] = 1.0;
    const auto result = nova::fem::estimate_diffusion_residual_error_3d(
        mesh, values, std::vector<double>{1.0, 1.0},
        [](double, double, double) { return 0.0; });

    const double face_term = 3.0 * std::sqrt(6.0) / 8.0;
    NOVA_REQUIRE(result.flux_jumps.faces.size() == 1U);
    NOVA_REQUIRE_NEAR(
        std::abs(result.flux_jumps.faces[0].normal_jump),
        std::sqrt(3.0) / 2.0, 2.0e-14);
    NOVA_REQUIRE_NEAR(result.element_errors[0], std::sqrt(face_term / 2.0), 2.0e-14);
    NOVA_REQUIRE_NEAR(result.element_errors[1], std::sqrt(face_term / 2.0), 2.0e-14);
    NOVA_REQUIRE_NEAR(result.total_error, std::sqrt(face_term), 2.0e-14);
}

NOVA_TEST(diffusion_residual_estimator_3d_detects_discontinuous_conductivity) {
    const auto mesh = two_tetrahedra_sharing_face();
    const H1FunctionSpace3D space(mesh);
    const auto values = space.interpolate(
        [](double x, double, double) { return x; });
    const auto result = nova::fem::estimate_diffusion_residual_error_3d(
        mesh, values, std::vector<double>{1.0, 2.0},
        [](double, double, double) { return 0.0; });

    const double face_term = std::sqrt(6.0) / 6.0;
    NOVA_REQUIRE(result.flux_jumps.faces.size() == 1U);
    NOVA_REQUIRE_NEAR(
        std::abs(result.flux_jumps.faces[0].normal_jump),
        1.0 / std::sqrt(3.0), 2.0e-14);
    NOVA_REQUIRE_NEAR(result.total_error, std::sqrt(face_term), 2.0e-14);
    NOVA_REQUIRE(result.element_errors[0] > 0.0);
    NOVA_REQUIRE(result.element_errors[1] > 0.0);
}

NOVA_TEST(diffusion_residual_estimator_3d_rejects_invalid_inputs) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const H1FunctionSpace3D space(box.mesh);
    auto values = space.interpolate(
        [](double x, double y, double z) { return x + y + z; });
    const std::vector<double> conductivity(box.mesh.elements.size(), 1.0);
    const auto zero_source = [](double, double, double) { return 0.0; };

    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::estimate_diffusion_residual_error_3d(
            box.mesh, std::vector<double>{0.0}, conductivity, zero_source);
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::estimate_diffusion_residual_error_3d(
            box.mesh, values, std::vector<double>{1.0}, zero_source);
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::estimate_diffusion_residual_error_3d(
            box.mesh, values,
            std::vector<double>(box.mesh.elements.size(), 0.0), zero_source);
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::estimate_diffusion_residual_error_3d(
            box.mesh, values, conductivity, {});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::estimate_diffusion_residual_error_3d(
            box.mesh, values, conductivity,
            [](double, double, double) {
                return std::numeric_limits<double>::infinity();
            });
    }));

    values[0] = std::numeric_limits<double>::quiet_NaN();
    NOVA_REQUIRE(rejects([&] {
        (void)nova::fem::estimate_diffusion_residual_error_3d(
            box.mesh, values, conductivity, zero_source);
    }));
}


// End-to-end regression across separate FEM solves. The estimator alone is not
// an error bound: compare both its global indicator and independently evaluated
// manufactured-solution errors after a geometrically conforming 1-to-8 refine.
NOVA_TEST(diffusion_residual_estimator_3d_solve_estimate_refine_resolve_converges) {
    using namespace nova::fem;

    const auto coarse_box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);

    const ExactScalarField3D exact{
        [](double x, double y, double z) { return x*x + y*y + z*z; },
        [](double x, double, double) { return 2.0*x; },
        [](double, double y, double) { return 2.0*y; },
        [](double, double, double z) { return 2.0*z; }
    };
    const auto source = [](double, double, double) { return -6.0; };

    ScalarDiffusion3DProblem coarse;
    coarse.mesh = coarse_box.mesh;
    coarse.cell_diffusivity.assign(coarse.mesh.elements.size(), 1.0);
    coarse.source = source;
    for (const auto& face : coarse_box.boundary_faces) {
        coarse.dirichlet_faces.push_back({face, exact.value});
    }

    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1.0e-12;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto coarse_solution = solve_scalar_diffusion_3d(coarse, options);
    NOVA_REQUIRE(coarse_solution.linear_solver.converged);

    const auto coarse_indicator = estimate_diffusion_residual_error_3d(
        coarse.mesh, coarse_solution.values, coarse.cell_diffusivity, source);
    const auto coarse_error = compute_error_norms_3d(
        coarse.mesh, coarse_solution.values, exact);
    NOVA_REQUIRE(std::isfinite(coarse_indicator.total_error));
    NOVA_REQUIRE(coarse_indicator.total_error > 0.0);
    NOVA_REQUIRE(coarse_error.l2 > 0.0);
    NOVA_REQUIRE(coarse_error.h1_seminorm > 0.0);

    const auto refinement = refine_tetrahedra_uniformly_3d(
        coarse.mesh, coarse_box.boundary_faces);
    NOVA_REQUIRE(refinement.mesh.elements.size() ==
                 8U * coarse.mesh.elements.size());
    NOVA_REQUIRE_NEAR(
        evaluate_mesh_quality_3d(refinement.mesh).total_volume,
        evaluate_mesh_quality_3d(coarse.mesh).total_volume, 1.0e-12);

    ScalarDiffusion3DProblem fine;
    fine.mesh = refinement.mesh;
    fine.cell_diffusivity = prolongate_cellwise_constant_uniform_3d(
        refinement, coarse.cell_diffusivity);
    fine.source = source;
    for (const auto& face : refinement.boundary_faces) {
        fine.dirichlet_faces.push_back({face, exact.value});
    }
    const auto fine_solution = solve_scalar_diffusion_3d(fine, options);
    NOVA_REQUIRE(fine_solution.linear_solver.converged);
    const auto fine_indicator = estimate_diffusion_residual_error_3d(
        fine.mesh, fine_solution.values, fine.cell_diffusivity, source);
    const auto fine_error = compute_error_norms_3d(
        fine.mesh, fine_solution.values, exact);

    NOVA_REQUIRE(fine_indicator.element_errors.size() == fine.mesh.elements.size());
    NOVA_REQUIRE(fine_indicator.total_error > 0.0);
    NOVA_REQUIRE(std::isfinite(fine_indicator.total_error));
    NOVA_REQUIRE(coarse_indicator.total_error > fine_indicator.total_error);
    NOVA_REQUIRE(coarse_error.l2 > fine_error.l2);
    NOVA_REQUIRE(coarse_error.h1_seminorm > fine_error.h1_seminorm);
}


NOVA_TEST(diffusion_residual_estimator_3d_dorfler_local_refine_and_resolve) {
    using namespace nova::fem;
    // Manufactured mixed-boundary Poisson problem: u=x(1-x), -Delta u=2.
    // Homogeneous Dirichlet x-end faces and natural zero-flux y/z faces
    // are physically consistent on EVERY refined mesh. The coarse P1
    // constrained trial space is nested in the refined trial space,
    // unlike re-interpolated nonhomogeneous quadratic boundary data.
    // Exercise estimator -> Dörfler -> LOCAL refine -> transfer -> re-solve.
    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const ExactScalarField3D exact{
        [](double x, double, double) { return x * (1.0 - x); },
        [](double x, double, double) { return 1.0 - 2.0*x; },
        [](double, double, double) { return 0.0; },
        [](double, double, double) { return 0.0; }
    };
    const auto source = [](double, double, double) { return 2.0; };
    const auto on_dirichlet_end = [](
        const TetrahedralMesh3D& mesh, const BoundaryFace3D& face
    ) {
        const double x = mesh.nodes[face[0]].x;
        return (x == 0.0 || x == 1.0) &&
            mesh.nodes[face[1]].x == x && mesh.nodes[face[2]].x == x;
    };

    ScalarDiffusion3DProblem coarse;
    coarse.mesh = box.mesh;
    coarse.cell_diffusivity.assign(coarse.mesh.elements.size(), 1.0);
    coarse.source = source;
    for (const auto& face : box.boundary_faces) {
        if (on_dirichlet_end(coarse.mesh, face)) {
            coarse.dirichlet_faces.push_back({face, exact.value});
        }
    }

    nova::linalg::SolverOptions options;
    options.relative_tolerance = 1.0e-12;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    const auto coarse_solution = solve_scalar_diffusion_3d(coarse, options);
    NOVA_REQUIRE(coarse_solution.linear_solver.converged);
    const auto before = compute_error_norms_3d(
        coarse.mesh, coarse_solution.values, exact);
    const auto indicators = estimate_diffusion_residual_error_3d(
        coarse.mesh, coarse_solution.values, coarse.cell_diffusivity, source);
    NOVA_REQUIRE(indicators.total_error > 0.0);

    const auto cycle = refine_tetrahedra_by_indicators_3d(
        coarse.mesh, indicators.element_errors, 0.65, box.boundary_faces);
    const auto& refinement = cycle.refinement;
    NOVA_REQUIRE(!cycle.marked_elements.empty());
    NOVA_REQUIRE(cycle.marked_elements.size() < coarse.mesh.elements.size());
    NOVA_REQUIRE(refinement.mesh.nodes.size() > coarse.mesh.nodes.size());
    NOVA_REQUIRE(refinement.mesh.elements.size() > coarse.mesh.elements.size());
    NOVA_REQUIRE(refinement.mesh.elements.size() <= 8U * coarse.mesh.elements.size());
    NOVA_REQUIRE(refinement.parent_elements.size() ==
                 refinement.mesh.elements.size());
    NOVA_REQUIRE(refinement.boundary_faces.size() ==
                 refinement.parent_boundary_faces.size());
    NOVA_REQUIRE_NEAR(
        evaluate_mesh_quality_3d(refinement.mesh).total_volume,
        evaluate_mesh_quality_3d(coarse.mesh).total_volume, 1.0e-12);

    // Transfer the previous solution and cell data without relying on the
    // solver to silently infer initial values or physical boundary metadata.
    const auto inherited = prolongate_p1_refined_3d(
        refinement, coarse_solution.values);
    NOVA_REQUIRE(inherited.size() == refinement.mesh.nodes.size());
    for (std::size_t i = 0; i < coarse.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(inherited[i], coarse_solution.values[i], 1.0e-14);
    }

    ScalarDiffusion3DProblem fine;
    fine.mesh = refinement.mesh;
    fine.cell_diffusivity = prolongate_cellwise_constant_refined_3d(
        refinement, coarse.cell_diffusivity);
    fine.source = source;
    for (const auto& face : refinement.boundary_faces) {
        if (on_dirichlet_end(fine.mesh, face)) {
            fine.dirichlet_faces.push_back({face, exact.value});
        }
    }

    const auto fine_solution = solve_scalar_diffusion_3d(fine, options);
    NOVA_REQUIRE(fine_solution.linear_solver.converged);
    const auto after = compute_error_norms_3d(
        fine.mesh, fine_solution.values, exact);
    NOVA_REQUIRE(std::isfinite(after.l2));
    NOVA_REQUIRE(std::isfinite(after.h1_seminorm));
    // The nested homogeneous-Dirichlet spaces make this a controlled
    // refinement test; no monotonicity claim is made for arbitrary physics.
    NOVA_REQUIRE(after.h1_seminorm < before.h1_seminorm);
    NOVA_REQUIRE(after.l2 < before.l2);
}


NOVA_TEST(diffusion_residual_estimator_3d_neumann_faces_match_outward_flux) {
    TetrahedralMesh3D mesh;
    mesh.nodes = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };
    mesh.elements = {{0U, 1U, 2U, 3U}};
    const H1FunctionSpace3D space(mesh);
    const auto values = space.interpolate(
        [](double x, double, double) { return x; });
    const std::vector<double> conductivity{2.0};
    const auto no_source = [](double, double, double) { return 0.0; };
    using nova::fem::NeumannResidualFace3D;
    std::vector<NeumannResidualFace3D> natural{
        {{0U, 2U, 3U}, [](double, double, double) { return 2.0; }}
    };

    const auto exact = nova::fem::estimate_diffusion_residual_error_3d(
        mesh, values, conductivity, no_source, natural);
    NOVA_REQUIRE(exact.boundary_residual_faces.size() == 1U);
    NOVA_REQUIRE(exact.boundary_residual_errors.size() == 1U);
    NOVA_REQUIRE_NEAR(exact.boundary_residual_faces[0].area, 0.5, 1.0e-14);
    NOVA_REQUIRE_NEAR(exact.boundary_residual_faces[0].diameter,
        std::sqrt(2.0), 1.0e-14);
    NOVA_REQUIRE(exact.total_error < 1.0e-12);

    // Input face winding must not alter the outward normal.
    natural[0].face = {3U, 2U, 0U};
    const auto reversed = nova::fem::estimate_diffusion_residual_error_3d(
        mesh, values, conductivity, no_source, natural);
    NOVA_REQUIRE(reversed.total_error < 1.0e-12);

    // g_N = 0 leaves a constant 2 W/m^2 flux mismatch on this face:
    // eta^2 = h_F * |F| * 2^2 = 2 * sqrt(2).
    natural[0].outward_flux =
        [](double, double, double) { return 0.0; };
    const auto nonzero = nova::fem::estimate_diffusion_residual_error_3d(
        mesh, values, conductivity, no_source, natural);
    const double expected = std::sqrt(2.0 * std::sqrt(2.0));
    NOVA_REQUIRE_NEAR(nonzero.boundary_residual_errors[0], expected, 3.0e-14);
    NOVA_REQUIRE_NEAR(nonzero.total_error, expected, 3.0e-14);
    NOVA_REQUIRE_NEAR(nonzero.boundary_residual_faces[0].mismatch_l2,
        std::sqrt(2.0), 3.0e-14);
}

NOVA_TEST(diffusion_residual_estimator_3d_integrates_varying_neumann_mismatch) {
    TetrahedralMesh3D mesh;
    mesh.nodes = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };
    mesh.elements = {{0U, 1U, 2U, 3U}};
    const std::vector<double> field(mesh.nodes.size(), 0.0);
    const std::vector<double> conductivity{1.0};
    const std::vector<nova::fem::NeumannResidualFace3D> natural{
        {{0U, 2U, 3U}, [](double, double y, double) { return y; }}
    };
    const auto result = nova::fem::estimate_diffusion_residual_error_3d(
        mesh, field, conductivity,
        [](double, double, double) { return 0.0; }, natural);

    // On x=0, integral_F y^2 dS = 1/12 exactly; h_F = sqrt(2).
    const double expected = std::sqrt(std::sqrt(2.0) / 12.0);
    NOVA_REQUIRE_NEAR(result.total_error, expected, 1.0e-13);
    NOVA_REQUIRE_NEAR(result.boundary_residual_errors[0], expected, 1.0e-13);
    NOVA_REQUIRE_NEAR(result.boundary_residual_faces[0].mismatch_l2,
        std::sqrt(1.0 / 12.0), 1.0e-13);
}

NOVA_TEST(diffusion_residual_estimator_3d_rejects_invalid_natural_faces) {
    const auto mesh = two_tetrahedra_sharing_face();
    const std::vector<double> field(mesh.nodes.size(), 0.0);
    const std::vector<double> conductivity(mesh.elements.size(), 1.0);
    const auto zero_source = [](double, double, double) { return 0.0; };
    using nova::fem::NeumannResidualFace3D;
    const NeumannResidualFace3D boundary{
        {0U, 1U, 2U}, [](double, double, double) { return 0.0; }};
    const auto check = [&](const std::vector<NeumannResidualFace3D>& natural) {
        return rejects([&] {
            (void)nova::fem::estimate_diffusion_residual_error_3d(
                mesh, field, conductivity, zero_source, natural);
        });
    };

    NOVA_REQUIRE(check({boundary, boundary}));
    NOVA_REQUIRE(check({{{1U, 2U, 3U},
        [](double, double, double) { return 0.0; }}})); // interior
    NOVA_REQUIRE(check({{{0U, 1U, mesh.nodes.size()},
        [](double, double, double) { return 0.0; }}})); // missing node
    NOVA_REQUIRE(check({{{0U, 1U, 4U},
        [](double, double, double) { return 0.0; }}})); // nonexistent face
    NOVA_REQUIRE(check({{{0U, 1U, 2U}, {}}})); // absent callback
    NOVA_REQUIRE(check({{{0U, 1U, 2U},
        [](double, double, double) {
            return std::numeric_limits<double>::quiet_NaN();
        }}}));

    // Backward-compatible path leaves out natural-boundary residuals.
    const auto old = nova::fem::estimate_diffusion_residual_error_3d(
        mesh, field, conductivity, zero_source);
    NOVA_REQUIRE(old.boundary_residual_faces.empty());
    NOVA_REQUIRE(old.boundary_residual_errors.size() == mesh.elements.size());
    NOVA_REQUIRE(old.total_error < 1.0e-12);
}


NOVA_TEST(diffusion_residual_estimator_3d_robin_matches_affine_field) {
    using namespace nova::fem;
    // Exact P1 field u=1+x, k=2, q=(-2,0,0).
    // On x=0, q.n=+2 (Neumann); on x=1, q.n=-2 =
    // h*(u-ambient)=2*(2-3) (Robin).
    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const auto values = H1FunctionSpace3D(box.mesh).interpolate(
        [](double x, double, double) { return 1.0 + x; });
    std::vector<NeumannResidualFace3D> neumann;
    std::vector<RobinResidualFace3D> robin;
    for (const auto& face : box.boundary_faces) {
        const double x = box.mesh.nodes[face[0]].x;
        if (box.mesh.nodes[face[1]].x != x ||
            box.mesh.nodes[face[2]].x != x) {
            continue;
        }
        if (x == 0.0) {
            neumann.push_back({face, [](double, double, double) { return 2.0; }});
        } else if (x == 1.0) {
            robin.push_back({
                face,
                [](double, double, double) { return 2.0; },
                [](double, double, double) { return 3.0; }
            });
        }
    }
    NOVA_REQUIRE(!neumann.empty());
    NOVA_REQUIRE(!robin.empty());
    const std::vector<double> conductivity(box.mesh.elements.size(), 2.0);
    const auto result = estimate_diffusion_residual_error_3d(
        box.mesh, values, conductivity,
        [](double, double, double) { return 0.0; },
        neumann, robin);
    NOVA_REQUIRE(result.total_error < 1.0e-12);
    NOVA_REQUIRE(result.boundary_residual_faces.size() ==
                 neumann.size() + robin.size());
    std::size_t robin_count = 0U;
    for (const auto& item : result.boundary_residual_faces) {
        NOVA_REQUIRE(item.mismatch_l2 < 1.0e-12);
        if (item.kind == BoundaryResidualKind3D::Robin) ++robin_count;
    }
    NOVA_REQUIRE(robin_count == robin.size());
}

NOVA_TEST(diffusion_residual_estimator_3d_robin_detects_wrong_ambient) {
    using namespace nova::fem;
    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const auto values = H1FunctionSpace3D(box.mesh).interpolate(
        [](double x, double, double) { return 1.0 + x; });
    std::vector<RobinResidualFace3D> robin;
    for (const auto& face : box.boundary_faces) {
        if (box.mesh.nodes[face[0]].x == 1.0 &&
            box.mesh.nodes[face[1]].x == 1.0 &&
            box.mesh.nodes[face[2]].x == 1.0) {
            robin.push_back({
                face,
                [](double, double, double) { return 2.0; },
                [](double, double, double) { return 0.0; }
            });
        }
    }
    const std::vector<double> conductivity(box.mesh.elements.size(), 2.0);
    const auto result = estimate_diffusion_residual_error_3d(
        box.mesh, values, conductivity,
        [](double, double, double) { return 0.0; },
        std::vector<NeumannResidualFace3D>{}, robin);
    // q.n=-2 but h*(u-ambient)=4 over the unit-area x=1 plate.
    NOVA_REQUIRE(result.total_error > 6.0);
    NOVA_REQUIRE(result.boundary_residual_faces.size() == robin.size());
    for (const auto& item : result.boundary_residual_faces) {
        NOVA_REQUIRE(item.kind == BoundaryResidualKind3D::Robin);
        NOVA_REQUIRE(item.mismatch_l2 > 0.0);
    }
}

NOVA_TEST(diffusion_residual_estimator_3d_robin_rejects_invalid_and_overlapping_faces) {
    using namespace nova::fem;
    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    const auto values = H1FunctionSpace3D(box.mesh).interpolate(
        [](double x, double, double) { return x; });
    const std::vector<double> conductivity(box.mesh.elements.size(), 2.0);
    const auto face = box.boundary_faces.front();
    const auto zero = [](double, double, double) { return 0.0; };
    const auto unit = [](double, double, double) { return 1.0; };
    const auto evaluate = [&](const std::vector<NeumannResidualFace3D>& neumann,
                              const std::vector<RobinResidualFace3D>& robin) {
        return estimate_diffusion_residual_error_3d(
            box.mesh, values, conductivity, zero, neumann, robin);
    };
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate({}, {{face, {}, zero}});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate({}, {{face, unit, {}}});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate({}, {{face,
            [](double, double, double) { return -1.0; }, zero}});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate({}, {{face, unit,
            [](double, double, double) {
                return std::numeric_limits<double>::infinity();
            }}});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate({}, {{face, unit, zero}, {face, unit, zero}});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)evaluate({{face, zero}}, {{face, unit, zero}});
    }));
}
