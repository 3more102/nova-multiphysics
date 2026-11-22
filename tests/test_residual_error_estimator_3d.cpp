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
    // Non-affine manufactured field: -Delta u = -6 with u = x^2+y^2+z^2.
    // Unlike the uniform-refinement regression, this uses the public
    // estimator -> Dörfler marking -> conforming LOCAL refinement -> P1/cell
    // transfer -> fresh FEM solve sequence on a genuinely nonzero residual.
    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    const ExactScalarField3D exact{
        [](double x, double y, double z) { return x*x + y*y + z*z; },
        [](double x, double, double) { return 2.0*x; },
        [](double, double y, double) { return 2.0*y; },
        [](double, double, double z) { return 2.0*z; }
    };
    const auto source = [](double, double, double) { return -6.0; };

    ScalarDiffusion3DProblem coarse;
    coarse.mesh = box.mesh;
    coarse.cell_diffusivity.assign(coarse.mesh.elements.size(), 1.0);
    coarse.source = source;
    for (const auto& face : box.boundary_faces) {
        coarse.dirichlet_faces.push_back({face, exact.value});
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
        fine.dirichlet_faces.push_back({face, exact.value});
    }

    const auto fine_solution = solve_scalar_diffusion_3d(fine, options);
    NOVA_REQUIRE(fine_solution.linear_solver.converged);
    const auto after = compute_error_norms_3d(
        fine.mesh, fine_solution.values, exact);
    NOVA_REQUIRE(std::isfinite(after.l2));
    NOVA_REQUIRE(std::isfinite(after.h1_seminorm));
    // These inequalities are test observations for this manufactured
    // problem, not a generic guarantee for arbitrary marked meshes.
    NOVA_REQUIRE(after.h1_seminorm < before.h1_seminorm);
    NOVA_REQUIRE(after.l2 < before.l2);
}
