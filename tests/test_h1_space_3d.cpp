#include "test_framework.hpp"

#include "nova/fem/h1_space_3d.hpp"
#include "nova/fem/tetrahedron_h1.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

NOVA_TEST(tetrahedron_p1_basis_has_partition_of_unity_and_zero_gradient_sum) {
    const auto basis =
        nova::fem::evaluate_tetrahedron_p1_basis_3d(0.2, 0.3, 0.1);
    double value_sum = 0.0;
    std::array<double, 3> gradient_sum{};
    for (std::size_t i = 0; i < 4U; ++i) {
        value_sum += basis.values[i];
        for (std::size_t d = 0; d < 3U; ++d) {
            gradient_sum[d] += basis.reference_gradients[i][d];
        }
    }
    NOVA_REQUIRE_NEAR(value_sum, 1.0, 1.0e-14);
    for (const double component : gradient_sum) {
        NOVA_REQUIRE_NEAR(component, 0.0, 1.0e-14);
    }
}

NOVA_TEST(h1_p1_3d_interpolates_affine_fields_and_recovers_physical_gradient) {
    const auto box = nova::fem::make_box_mesh_3d(
        -1.0, 2.0, 0.0, 2.0, -0.5, 1.5, 2U, 2U, 2U
    );
    const nova::fem::H1FunctionSpace3D space(box.mesh);
    const auto field = [](const double x, const double y, const double z) {
        return 1.5 + 2.0 * x - 3.0 * y + 4.0 * z;
    };
    const auto values = space.interpolate(field);
    NOVA_REQUIRE(values.size() == space.dimension());

    for (std::size_t element = 0; element < box.mesh.elements.size(); ++element) {
        const auto physical =
            nova::fem::evaluate_tetrahedron_p1_physical_gradients_3d(
                box.mesh, element
            );
        std::array<double, 3> gradient{};
        const auto& dofs = space.dofs().element_dofs(element);
        for (std::size_t local = 0; local < 4U; ++local) {
            for (std::size_t d = 0; d < 3U; ++d) {
                gradient[d] += values[dofs[local]] *
                    physical.gradients[local][d];
            }
        }
        NOVA_REQUIRE_NEAR(gradient[0], 2.0, 2.0e-12);
        NOVA_REQUIRE_NEAR(gradient[1], -3.0, 2.0e-12);
        NOVA_REQUIRE_NEAR(gradient[2], 4.0, 2.0e-12);
    }
}

NOVA_TEST(h1_p1_3d_handles_reversed_tetrahedron_orientation) {
    auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U
    );
    std::swap(box.mesh.elements[0][1], box.mesh.elements[0][3]);
    const nova::fem::H1FunctionSpace3D space(box.mesh);
    const auto values = space.interpolate(
        [](double x, double y, double z) { return x + 2.0 * y - z; }
    );
    const auto physical =
        nova::fem::evaluate_tetrahedron_p1_physical_gradients_3d(
            box.mesh, 0U
        );
    std::array<double, 3> gradient{};
    const auto& dofs = space.dofs().element_dofs(0U);
    for (std::size_t local = 0; local < 4U; ++local) {
        for (std::size_t d = 0; d < 3U; ++d) {
            gradient[d] += values[dofs[local]] * physical.gradients[local][d];
        }
    }
    NOVA_REQUIRE_NEAR(gradient[0], 1.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(gradient[1], 2.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(gradient[2], -1.0, 1.0e-12);
}

NOVA_TEST(h1_p1_3d_interpolation_rejects_missing_or_nonfinite_field) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U
    );
    const nova::fem::H1FunctionSpace3D space(box.mesh);

    bool rejected_empty = false;
    try {
        static_cast<void>(space.interpolate({}));
    } catch (const std::invalid_argument&) {
        rejected_empty = true;
    }
    NOVA_REQUIRE(rejected_empty);

    bool rejected_nonfinite = false;
    try {
        static_cast<void>(space.interpolate(
            [](double, double, double) {
                return std::numeric_limits<double>::quiet_NaN();
            }
        ));
    } catch (const std::invalid_argument&) {
        rejected_nonfinite = true;
    }
    NOVA_REQUIRE(rejected_nonfinite);
}
