#include "test_framework.hpp"

#include "nova/fem/error_estimator_2d.hpp"
#include "nova/fem/h1_space_2d.hpp"
#include "nova/fem/mesh_refinement_2d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

NOVA_TEST(mesh_uniform_refinement_prolongates_quadratic_p2_exactly) {
    auto rectangle = nova::fem::make_rectangle_mesh_2d(
        -0.5, 1.25, -0.75, 1.5, 3U, 2U
    );
    std::swap(rectangle.mesh.elements[0][1], rectangle.mesh.elements[0][2]);

    const auto quadratic = [](const double x, const double y) {
        return 1.25 + 0.7 * x - 1.1 * y +
            0.3 * x * x + 0.4 * x * y - 0.2 * y * y;
    };

    const nova::fem::H1FunctionSpace2D coarse_space(
        rectangle.mesh, nova::fem::H1TriangleOrder::P2
    );
    const auto coarse_values = coarse_space.interpolate(quadratic);
    const auto refined = nova::fem::refine_triangles_uniformly_2d(
        rectangle.mesh, rectangle.boundary_edges
    );
    const auto fine_values = nova::fem::prolongate_p2_uniform_2d(
        rectangle.mesh, refined, coarse_values
    );
    const nova::fem::H1FunctionSpace2D fine_space(
        refined.mesh, nova::fem::H1TriangleOrder::P2
    );
    const auto exact_fine_values = fine_space.interpolate(quadratic);

    NOVA_REQUIRE(fine_values.size() == fine_space.dimension());
    NOVA_REQUIRE(fine_values.size() == exact_fine_values.size());
    for (std::size_t dof = 0; dof < fine_values.size(); ++dof) {
        NOVA_REQUIRE_NEAR(fine_values[dof], exact_fine_values[dof], 2.0e-12);
    }
}

NOVA_TEST(mesh_uniform_refinement_p2_transfer_rejects_invalid_metadata_and_state) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1U, 1U
    );
    const nova::fem::H1FunctionSpace2D coarse_space(
        rectangle.mesh, nova::fem::H1TriangleOrder::P2
    );
    const auto field = coarse_space.interpolate(
        [](const double x, const double y) { return x + y; }
    );
    const auto refined = nova::fem::refine_triangles_uniformly_2d(rectangle.mesh);

    bool rejected_size = false;
    try {
        static_cast<void>(nova::fem::prolongate_p2_uniform_2d(
            rectangle.mesh, refined, std::vector<double>(1U, 0.0)
        ));
    } catch (const std::invalid_argument&) {
        rejected_size = true;
    }
    NOVA_REQUIRE(rejected_size);

    auto non_finite = field;
    non_finite[0] = std::numeric_limits<double>::quiet_NaN();
    bool rejected_non_finite = false;
    try {
        static_cast<void>(nova::fem::prolongate_p2_uniform_2d(
            rectangle.mesh, refined, non_finite
        ));
    } catch (const std::invalid_argument&) {
        rejected_non_finite = true;
    }
    NOVA_REQUIRE(rejected_non_finite);

    auto invalid_parent = refined;
    invalid_parent.parent_elements[0] = rectangle.mesh.elements.size();
    bool rejected_parent = false;
    try {
        static_cast<void>(nova::fem::prolongate_p2_uniform_2d(
            rectangle.mesh, invalid_parent, field
        ));
    } catch (const std::out_of_range&) {
        rejected_parent = true;
    }
    NOVA_REQUIRE(rejected_parent);
}

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
    const auto estimate = nova::fem::estimate_p1_diffusion_error_2d(
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
    const auto estimate = nova::fem::estimate_p1_diffusion_error_2d(
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
    const auto estimate = nova::fem::estimate_p1_diffusion_error_2d(
        rectangle.mesh, values, diffusivity,
        [](const double, const double) { return 0.0; }
    );

    NOVA_REQUIRE(estimate.element_indicators.size() == 2U);
    NOVA_REQUIRE_NEAR(estimate.element_indicators[0], std::sqrt(0.5), 1.0e-13);
    NOVA_REQUIRE_NEAR(estimate.element_indicators[1], std::sqrt(0.5), 1.0e-13);
    NOVA_REQUIRE_NEAR(estimate.global_estimate, 1.0, 1.0e-13);
}

NOVA_TEST(dorfler_marking_selects_largest_bulk_indicators_deterministically) {
    const std::vector<double> indicators{1.0, 2.0, 3.0, 4.0};
    const auto half = nova::fem::mark_elements_dorfler_2d(indicators, 0.5);
    NOVA_REQUIRE(half.size() == 1U);
    NOVA_REQUIRE(half[0] == 3U);

    const auto eighty_percent =
        nova::fem::mark_elements_dorfler_2d(indicators, 0.8);
    NOVA_REQUIRE(eighty_percent.size() == 2U);
    NOVA_REQUIRE(eighty_percent[0] == 3U);
    NOVA_REQUIRE(eighty_percent[1] == 2U);

    const std::vector<double> zero{0.0, 0.0, 0.0};
    NOVA_REQUIRE(nova::fem::mark_elements_dorfler_2d(zero, 0.5).empty());
}

NOVA_TEST(p1_residual_estimator_and_dorfler_marking_reject_invalid_inputs) {
    const nova::fem::TriangleMesh2D mesh{{
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}
    }, {{0U, 1U, 2U}}};
    const std::vector<double> values{0.0, 0.0, 0.0};
    const std::vector<double> diffusivity{1.0};

    bool rejected_diffusivity = false;
    try {
        const std::vector<double> invalid_k{0.0};
        static_cast<void>(nova::fem::estimate_p1_diffusion_error_2d(
            mesh, values, invalid_k,
            [](const double, const double) { return 0.0; }
        ));
    } catch (const std::invalid_argument&) {
        rejected_diffusivity = true;
    }
    NOVA_REQUIRE(rejected_diffusivity);

    bool rejected_source = false;
    try {
        static_cast<void>(nova::fem::estimate_p1_diffusion_error_2d(
            mesh, values, diffusivity,
            [](const double, const double) {
                return std::numeric_limits<double>::infinity();
            }
        ));
    } catch (const std::invalid_argument&) {
        rejected_source = true;
    }
    NOVA_REQUIRE(rejected_source);

    bool rejected_theta = false;
    try {
        static_cast<void>(nova::fem::mark_elements_dorfler_2d(
            std::vector<double>{1.0}, 0.0
        ));
    } catch (const std::invalid_argument&) {
        rejected_theta = true;
    }
    NOVA_REQUIRE(rejected_theta);
}
