#include "test_framework.hpp"

#include "nova/fem/h1_space_2d.hpp"

#include <cmath>
#include <stdexcept>

namespace {

nova::fem::TriangleMesh2D make_two_triangle_mesh() {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0}
    };
    mesh.elements = {
        {0, 1, 2},
        {0, 2, 3}
    };
    return mesh;
}

void require_partition_of_unity(
    const nova::fem::H1TriangleOrder order,
    const double r,
    const double s
) {
    const auto basis = nova::fem::evaluate_triangle_h1_basis(order, r, s);

    double value_sum = 0.0;
    double gradient_x_sum = 0.0;
    double gradient_y_sum = 0.0;

    for (std::size_t i = 0; i < basis.count; ++i) {
        value_sum += basis.values[i];
        gradient_x_sum += basis.reference_gradients[i][0];
        gradient_y_sum += basis.reference_gradients[i][1];
    }

    NOVA_REQUIRE_NEAR(value_sum, 1.0, 1.0e-13);
    NOVA_REQUIRE_NEAR(gradient_x_sum, 0.0, 1.0e-13);
    NOVA_REQUIRE_NEAR(gradient_y_sum, 0.0, 1.0e-13);
}

} // namespace

NOVA_TEST(h1_triangle_p1_basis_has_partition_of_unity) {
    require_partition_of_unity(
        nova::fem::H1TriangleOrder::P1,
        0.23,
        0.31
    );
}

NOVA_TEST(h1_triangle_p2_basis_has_partition_of_unity) {
    require_partition_of_unity(
        nova::fem::H1TriangleOrder::P2,
        0.23,
        0.31
    );
}

NOVA_TEST(h1_triangle_p2_basis_has_nodal_kronecker_property) {
    using nova::fem::H1TriangleOrder;
    using nova::fem::evaluate_triangle_h1_basis;

    const auto v0 = evaluate_triangle_h1_basis(H1TriangleOrder::P2, 0.0, 0.0);
    const auto v1 = evaluate_triangle_h1_basis(H1TriangleOrder::P2, 1.0, 0.0);
    const auto v2 = evaluate_triangle_h1_basis(H1TriangleOrder::P2, 0.0, 1.0);
    const auto e01 = evaluate_triangle_h1_basis(H1TriangleOrder::P2, 0.5, 0.0);
    const auto e12 = evaluate_triangle_h1_basis(H1TriangleOrder::P2, 0.5, 0.5);
    const auto e20 = evaluate_triangle_h1_basis(H1TriangleOrder::P2, 0.0, 0.5);

    const std::array<std::array<double, 6>, 6> nodal_values{{
        v0.values, v1.values, v2.values, e01.values, e12.values, e20.values
    }};

    for (std::size_t row = 0; row < nodal_values.size(); ++row) {
        for (std::size_t column = 0; column < nodal_values[row].size(); ++column) {
            const double expected = row == column ? 1.0 : 0.0;
            NOVA_REQUIRE_NEAR(nodal_values[row][column], expected, 1.0e-13);
        }
    }
}

NOVA_TEST(h1_dof_manager_p1_uses_vertex_dofs_only) {
    const auto mesh = make_two_triangle_mesh();
    const nova::fem::H1DofManager2D dofs(
        mesh,
        nova::fem::H1TriangleOrder::P1
    );

    NOVA_REQUIRE(dofs.dof_count() == 4);
    NOVA_REQUIRE(dofs.vertex_dof_count() == 4);
    NOVA_REQUIRE(dofs.edge_dof_count() == 5);
    NOVA_REQUIRE(dofs.element_dof_count() == 3);
    NOVA_REQUIRE(dofs.element_dofs(0)[0] == 0);
    NOVA_REQUIRE(dofs.element_dofs(0)[1] == 1);
    NOVA_REQUIRE(dofs.element_dofs(0)[2] == 2);
    NOVA_REQUIRE(!dofs.edge_dof(0, 1).has_value());
}

NOVA_TEST(h1_dof_manager_p2_shares_edge_dofs_between_elements) {
    const auto mesh = make_two_triangle_mesh();
    const nova::fem::H1DofManager2D dofs(
        mesh,
        nova::fem::H1TriangleOrder::P2
    );

    NOVA_REQUIRE(dofs.dof_count() == 9);
    NOVA_REQUIRE(dofs.vertex_dof_count() == 4);
    NOVA_REQUIRE(dofs.edge_dof_count() == 5);
    NOVA_REQUIRE(dofs.element_dof_count() == 6);

    const auto diagonal = dofs.edge_dof(0, 2);
    const auto reversed_diagonal = dofs.edge_dof(2, 0);
    NOVA_REQUIRE(diagonal.has_value());
    NOVA_REQUIRE(reversed_diagonal.has_value());
    NOVA_REQUIRE(diagonal.value() == reversed_diagonal.value());

    NOVA_REQUIRE(dofs.element_dofs(0)[5] == diagonal.value());
    NOVA_REQUIRE(dofs.element_dofs(1)[3] == diagonal.value());
    NOVA_REQUIRE(dofs.dof_entity(diagonal.value()) == nova::fem::H1DofEntity::Edge);
    NOVA_REQUIRE(dofs.dof_entity(0) == nova::fem::H1DofEntity::Vertex);
}

NOVA_TEST(h1_function_space_p1_interpolates_vertices) {
    const auto mesh = make_two_triangle_mesh();
    const nova::fem::H1FunctionSpace2D space(
        mesh,
        nova::fem::H1TriangleOrder::P1
    );

    const auto values = space.interpolate(
        [](const double x, const double y) { return 2.0 * x - y + 3.0; }
    );

    NOVA_REQUIRE(values.size() == 4);
    NOVA_REQUIRE_NEAR(values[0], 3.0, 1.0e-13);
    NOVA_REQUIRE_NEAR(values[1], 5.0, 1.0e-13);
    NOVA_REQUIRE_NEAR(values[2], 4.0, 1.0e-13);
    NOVA_REQUIRE_NEAR(values[3], 2.0, 1.0e-13);
}

NOVA_TEST(h1_function_space_p2_interpolates_edge_midpoints) {
    const auto mesh = make_two_triangle_mesh();
    const nova::fem::H1FunctionSpace2D space(
        mesh,
        nova::fem::H1TriangleOrder::P2
    );

    const auto values = space.interpolate(
        [](const double x, const double y) { return x * x + 2.0 * y; }
    );

    NOVA_REQUIRE(values.size() == 9);

    for (std::size_t edge = 0; edge < space.dofs().edge_keys().size(); ++edge) {
        const auto [first, second] = space.dofs().edge_keys()[edge];
        const auto& p0 = mesh.nodes[first];
        const auto& p1 = mesh.nodes[second];
        const double x = 0.5 * (p0.x + p1.x);
        const double y = 0.5 * (p0.y + p1.y);
        NOVA_REQUIRE_NEAR(
            values[mesh.nodes.size() + edge],
            x * x + 2.0 * y,
            1.0e-13
        );
    }
}

NOVA_TEST(h1_dof_manager_rejects_invalid_element_index) {
    const auto mesh = make_two_triangle_mesh();
    const nova::fem::H1DofManager2D dofs(
        mesh,
        nova::fem::H1TriangleOrder::P2
    );

    bool threw = false;
    try {
        static_cast<void>(dofs.element_dofs(2));
    } catch (const std::out_of_range&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
