#include "test_framework.hpp"

#include "nova/fem/error_estimator_2d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

NOVA_TEST(p1_residual_estimator_is_zero_for_affine_harmonic_field) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 2.0, -1.0, 1.0, 3U, 2U
    );
    std::vector<double> values;
    for (const auto& node : rectangle.mesh.nodes) {
        values.push_back(1.0 + 2.0 * node.x - 3.0 * node.y);
    }
    const std::vector<double> diffusivity(
        rectangle.mesh.elements.size(), 2.5
    );
    const auto estimate = nova::fem::estimate_p1_residual_indicators_2d(
        rectangle.mesh, values, diffusivity,
        [](const double, const double) { return 0.0; }
    );

    NOVA_REQUIRE(estimate.element_indicators.size() ==
        rectangle.mesh.elements.size());
    NOVA_REQUIRE_NEAR(estimate.global_estimate, 0.0, 1.0e-13);
}

NOVA_TEST(p1_residual_estimator_integrates_constant_cell_residual) {
    const nova::fem::TriangleMesh2D mesh{{
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}
    }, {{0U, 1U, 2U}}};
    const std::vector<double> values{0.0, 0.0, 0.0};
    const std::vector<double> diffusivity{1.0};
    const auto estimate = nova::fem::estimate_p1_residual_indicators_2d(
        mesh, values, diffusivity,
        [](const double, const double) { return 3.0; }
    );

    NOVA_REQUIRE(estimate.element_indicators.size() == 1U);
    NOVA_REQUIRE_NEAR(estimate.element_indicators[0], 3.0, 1.0e-13);
    NOVA_REQUIRE_NEAR(estimate.global_estimate, 3.0, 1.0e-13);
}

NOVA_TEST(p1_residual_estimator_detects_interior_flux_jump) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    const std::vector<double> values{0.0, 1.0, 0.0, 1.0};
    const std::vector<double> diffusivity{1.0, 2.0};
    const auto estimate = nova::fem::estimate_p1_residual_indicators_2d(
        rectangle.mesh, values, diffusivity,
        [](const double, const double) { return 0.0; }
    );

    NOVA_REQUIRE(estimate.element_indicators.size() == 2U);
    NOVA_REQUIRE_NEAR(estimate.element_indicators[0], std::sqrt(0.5), 1.0e-13);
    NOVA_REQUIRE_NEAR(estimate.element_indicators[1], std::sqrt(0.5), 1.0e-13);
    NOVA_REQUIRE_NEAR(estimate.global_estimate, 1.0, 1.0e-13);
}

NOVA_TEST(p1_residual_estimator_rejects_invalid_inputs) {
    const nova::fem::TriangleMesh2D mesh{{
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}
    }, {{0U, 1U, 2U}}};
    const std::vector<double> values{0.0, 0.0, 0.0};
    const std::vector<double> diffusivity{1.0};

    bool rejected_size = false;
    try {
        static_cast<void>(nova::fem::estimate_p1_residual_indicators_2d(
            mesh, std::vector<double>{0.0}, diffusivity,
            [](const double, const double) { return 0.0; }
        ));
    } catch (const std::invalid_argument&) {
        rejected_size = true;
    }
    NOVA_REQUIRE(rejected_size);

    bool rejected_diffusivity = false;
    try {
        const std::vector<double> invalid_k{0.0};
        static_cast<void>(nova::fem::estimate_p1_residual_indicators_2d(
            mesh, values, invalid_k,
            [](const double, const double) { return 0.0; }
        ));
    } catch (const std::invalid_argument&) {
        rejected_diffusivity = true;
    }
    NOVA_REQUIRE(rejected_diffusivity);

    bool rejected_source = false;
    try {
        static_cast<void>(nova::fem::estimate_p1_residual_indicators_2d(
            mesh, values, diffusivity,
            [](const double, const double) {
                return std::numeric_limits<double>::infinity();
            }
        ));
    } catch (const std::invalid_argument&) {
        rejected_source = true;
    }
    NOVA_REQUIRE(rejected_source);
}
