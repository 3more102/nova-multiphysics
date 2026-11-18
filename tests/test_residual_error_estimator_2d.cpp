#include "test_framework.hpp"

#include "nova/fem/h1_space_2d.hpp"
#include "nova/fem/residual_error_estimator_2d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using nova::fem::H1FunctionSpace2D;
using nova::fem::H1TriangleOrder;
using nova::fem::TriangleMesh2D;

} // namespace

NOVA_TEST(diffusion_residual_estimator_is_zero_for_affine_harmonic_p1_field) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 2U, 2U
    );
    const H1FunctionSpace2D space(rectangle.mesh, H1TriangleOrder::P1);
    const auto values = space.interpolate([](const double x, const double y) {
        return 1.0 + 2.0 * x - 3.0 * y;
    });
    const std::vector<double> conductivity(rectangle.mesh.elements.size(), 1.0);
    const auto result = nova::fem::estimate_diffusion_residual_error_2d(
        rectangle.mesh, H1TriangleOrder::P1, values, conductivity,
        [](const double, const double) { return 0.0; }
    );
    NOVA_REQUIRE(result.element_errors.size() == rectangle.mesh.elements.size());
    NOVA_REQUIRE(result.total_error < 1.0e-12);
}

NOVA_TEST(diffusion_residual_estimator_cancels_exact_p2_cell_residual) {
    TriangleMesh2D mesh;
    mesh.nodes = {{0.0, 0.0}, {2.0, 0.0}, {1.0, 1.0}};
    mesh.elements = {{0U, 1U, 2U}};
    const H1FunctionSpace2D space(mesh, H1TriangleOrder::P2);
    const auto values = space.interpolate([](const double x, const double y) {
        return x * x + y * y;
    });
    const std::vector<double> conductivity{1.0};
    const auto result = nova::fem::estimate_diffusion_residual_error_2d(
        mesh, H1TriangleOrder::P2, values, conductivity,
        [](const double, const double) { return -4.0; }
    );
    NOVA_REQUIRE(result.element_errors.size() == 1U);
    NOVA_REQUIRE(result.total_error < 2.0e-12);
}

NOVA_TEST(diffusion_residual_estimator_integrates_constant_cell_source) {
    TriangleMesh2D mesh;
    mesh.nodes = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    mesh.elements = {{0U, 1U, 2U}};
    const H1FunctionSpace2D space(mesh, H1TriangleOrder::P1);
    const auto values = space.interpolate(
        [](const double, const double) { return 0.0; }
    );
    const auto result = nova::fem::estimate_diffusion_residual_error_2d(
        mesh, H1TriangleOrder::P1, values, std::vector<double>{1.0},
        [](const double, const double) { return 1.0; }
    );
    NOVA_REQUIRE_NEAR(result.element_errors[0], 1.0, 1.0e-13);
    NOVA_REQUIRE_NEAR(result.total_error, 1.0, 1.0e-13);
}

NOVA_TEST(diffusion_residual_estimator_detects_discontinuous_conductivity_flux_jump) {
    TriangleMesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}
    };
    mesh.elements = {{0U, 1U, 2U}, {0U, 2U, 3U}};
    const H1FunctionSpace2D space(mesh, H1TriangleOrder::P1);
    const auto values = space.interpolate(
        [](const double x, const double) { return x; }
    );
    const auto zero_source = [](const double, const double) { return 0.0; };

    const auto continuous_flux = nova::fem::estimate_diffusion_residual_error_2d(
        mesh, H1TriangleOrder::P1, values,
        std::vector<double>{1.0, 1.0}, zero_source
    );
    NOVA_REQUIRE(continuous_flux.total_error < 1.0e-12);

    // Reverse one cell to verify that the outward-normal convention is
    // independent of mixed clockwise/counter-clockwise element orientation.
    std::swap(mesh.elements[1][1], mesh.elements[1][2]);
    const H1FunctionSpace2D mixed_space(mesh, H1TriangleOrder::P1);
    const auto mixed_values = mixed_space.interpolate(
        [](const double x, const double) { return x; }
    );
    const auto mixed_continuous = nova::fem::estimate_diffusion_residual_error_2d(
        mesh, H1TriangleOrder::P1, mixed_values,
        std::vector<double>{1.0, 1.0}, zero_source
    );
    NOVA_REQUIRE(mixed_continuous.total_error < 1.0e-12);

    const auto discontinuous_flux =
        nova::fem::estimate_diffusion_residual_error_2d(
            mesh, H1TriangleOrder::P1, mixed_values,
            std::vector<double>{1.0, 2.0}, zero_source
        );
    NOVA_REQUIRE(discontinuous_flux.total_error > 0.1);
    NOVA_REQUIRE(discontinuous_flux.element_errors[0] > 0.0);
    NOVA_REQUIRE(discontinuous_flux.element_errors[1] > 0.0);
}

NOVA_TEST(diffusion_residual_estimator_rejects_invalid_inputs) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    const H1FunctionSpace2D space(rectangle.mesh, H1TriangleOrder::P2);
    auto values = space.interpolate(
        [](const double x, const double y) { return x + y; }
    );
    const std::vector<double> conductivity(rectangle.mesh.elements.size(), 1.0);
    const auto zero_source = [](const double, const double) { return 0.0; };

    bool rejected_size = false;
    try {
        static_cast<void>(nova::fem::estimate_diffusion_residual_error_2d(
            rectangle.mesh, H1TriangleOrder::P2,
            std::vector<double>{0.0}, conductivity, zero_source
        ));
    } catch (const std::invalid_argument&) {
        rejected_size = true;
    }
    NOVA_REQUIRE(rejected_size);

    bool rejected_conductivity_size = false;
    try {
        static_cast<void>(nova::fem::estimate_diffusion_residual_error_2d(
            rectangle.mesh, H1TriangleOrder::P2, values,
            std::vector<double>{1.0}, zero_source
        ));
    } catch (const std::invalid_argument&) {
        rejected_conductivity_size = true;
    }
    NOVA_REQUIRE(rejected_conductivity_size);

    bool rejected_nonpositive_conductivity = false;
    try {
        static_cast<void>(nova::fem::estimate_diffusion_residual_error_2d(
            rectangle.mesh, H1TriangleOrder::P2, values,
            std::vector<double>(rectangle.mesh.elements.size(), 0.0),
            zero_source
        ));
    } catch (const std::invalid_argument&) {
        rejected_nonpositive_conductivity = true;
    }
    NOVA_REQUIRE(rejected_nonpositive_conductivity);

    bool rejected_missing_source = false;
    try {
        static_cast<void>(nova::fem::estimate_diffusion_residual_error_2d(
            rectangle.mesh, H1TriangleOrder::P2, values, conductivity, {}
        ));
    } catch (const std::invalid_argument&) {
        rejected_missing_source = true;
    }
    NOVA_REQUIRE(rejected_missing_source);

    values[0] = std::numeric_limits<double>::quiet_NaN();
    bool rejected_nonfinite_state = false;
    try {
        static_cast<void>(nova::fem::estimate_diffusion_residual_error_2d(
            rectangle.mesh, H1TriangleOrder::P2, values, conductivity, zero_source
        ));
    } catch (const std::invalid_argument&) {
        rejected_nonfinite_state = true;
    }
    NOVA_REQUIRE(rejected_nonfinite_state);
}
