#include "test_framework.hpp"

#include "nova/fem/error_indicators_2d.hpp"
#include "nova/fem/h1_space_2d.hpp"
#include "nova/fem/mesh_transfer_2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

} // namespace

NOVA_TEST(gradient_recovery_is_exact_for_global_quadratic_p2_field) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        -1.0, 1.0, -0.5, 1.5, 3U, 2U
    );
    const auto field = [](const double x, const double y) {
        return 1.0 + 2.0 * x - y + x * x + 0.5 * x * y + 2.0 * y * y;
    };
    const nova::fem::H1FunctionSpace2D space(
        rectangle.mesh, nova::fem::H1TriangleOrder::P2
    );
    const auto values = space.interpolate(field);
    const auto indicators =
        nova::fem::estimate_h1_gradient_recovery_error_2d(
            rectangle.mesh, nova::fem::H1TriangleOrder::P2, values
        );

    NOVA_REQUIRE(indicators.element_errors.size() ==
        rectangle.mesh.elements.size());
    NOVA_REQUIRE(indicators.total_error < 2.0e-11);
}

NOVA_TEST(gradient_recovery_dorfler_marking_selects_largest_error_cells) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 3U, 3U
    );
    const nova::fem::H1FunctionSpace2D space(
        rectangle.mesh, nova::fem::H1TriangleOrder::P1
    );
    const auto values = space.interpolate(
        [](const double x, const double y) { return x * y; }
    );
    const auto indicators =
        nova::fem::estimate_h1_gradient_recovery_error_2d(
            rectangle.mesh, nova::fem::H1TriangleOrder::P1, values
        );
    NOVA_REQUIRE(indicators.total_error > 0.0);

    const double fraction = 0.6;
    const auto marked = nova::fem::mark_elements_dorfler_2d(
        indicators.element_errors, fraction
    );
    NOVA_REQUIRE(!marked.empty());
    NOVA_REQUIRE(std::is_sorted(marked.begin(), marked.end()));

    double selected_squared = 0.0;
    for (const std::size_t element : marked) {
        const double eta = indicators.element_errors[element];
        selected_squared += eta * eta;
    }
    NOVA_REQUIRE(
        selected_squared + 1.0e-14 >=
        fraction * indicators.total_error * indicators.total_error
    );
}

NOVA_TEST(transient_heat_gradient_recovery_uses_problem_conductivity) {
    const auto rectangle = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 2U, 2U
    );
    nova::fem::TransientHeat2DProblem problem;
    problem.mesh = rectangle.mesh;
    problem.conductivity = [](const double x, const double y) {
        return 1.0 + x + y;
    };
    const nova::fem::H1FunctionSpace2D space(
        problem.mesh, nova::fem::H1TriangleOrder::P1
    );
    const auto values = space.interpolate(
        [](const double x, const double y) { return x * y; }
    );
    const auto from_problem =
        nova::fem::estimate_transient_heat_gradient_recovery_2d(
            problem, nova::fem::H1TriangleOrder::P1, values
        );
    const auto direct =
        nova::fem::estimate_h1_gradient_recovery_error_2d(
            problem.mesh, nova::fem::H1TriangleOrder::P1, values,
            problem.conductivity
        );
    NOVA_REQUIRE_NEAR(
        from_problem.total_error, direct.total_error, 1.0e-14
    );
}

NOVA_TEST(dorfler_marking_rejects_invalid_inputs_and_handles_zero_error) {
    const std::array<double, 3> zeros{{0.0, 0.0, 0.0}};
    NOVA_REQUIRE(
        nova::fem::mark_elements_dorfler_2d(zeros, 0.5).empty()
    );

    bool rejected_fraction = false;
    try {
        static_cast<void>(nova::fem::mark_elements_dorfler_2d(zeros, 0.0));
    } catch (const std::invalid_argument&) {
        rejected_fraction = true;
    }
    NOVA_REQUIRE(rejected_fraction);

    const std::array<double, 2> negative{{1.0, -0.1}};
    bool rejected_error = false;
    try {
        static_cast<void>(
            nova::fem::mark_elements_dorfler_2d(negative, 0.5)
        );
    } catch (const std::invalid_argument&) {
        rejected_error = true;
    }
    NOVA_REQUIRE(rejected_error);
}
