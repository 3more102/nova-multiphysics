#include "test_framework.hpp"

#include "nova/fem/diffusion_2d.hpp"
#include "nova/fem/mesh2d.hpp"

#include <stdexcept>
#include <utility>

NOVA_TEST(mesh2d_rejects_non_manifold_edges) {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
        {0.0, -1.0},
        {0.5, 0.5}
    };
    mesh.elements = {
        {0, 1, 2},
        {1, 0, 3},
        {0, 1, 4}
    };

    bool threw = false;
    try {
        mesh.validate();
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(mesh2d_rejects_duplicate_triangle_with_reversed_connectivity) {
    nova::fem::TriangleMesh2D mesh;
    mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0}
    };
    mesh.elements = {
        {0, 1, 2},
        {2, 1, 0}
    };

    bool threw = false;
    try {
        mesh.validate();
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}

NOVA_TEST(diffusion_2d_rejects_duplicate_natural_edge_conditions) {
    auto generated = nova::fem::make_rectangle_mesh_2d(
        0.0, 1.0, 0.0, 1.0, 1, 1
    );

    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(generated.mesh);
    problem.diffusivity = [](const double, const double) { return 1.0; };
    problem.source = [](const double, const double) { return 0.0; };

    for (const std::size_t node : generated.boundary_nodes) {
        problem.dirichlet_constraints.push_back({node, 0.0});
    }

    problem.neumann_conditions.push_back({
        {0, 1},
        [](const double, const double) { return 0.0; }
    });
    problem.robin_conditions.push_back({
        {1, 0},
        [](const double, const double) { return 1.0; },
        [](const double, const double) { return 0.0; }
    });

    bool threw = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_2d(problem));
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    NOVA_REQUIRE(threw);
}
