#include "test_framework.hpp"

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
