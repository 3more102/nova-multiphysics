#include "test_framework.hpp"

#include "nova/fem/diffusion_2d.hpp"
#include "nova/fem/h1_space_2d.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace {

void constrain_vertical_boundaries(
    nova::fem::ScalarDiffusion2DProblem& problem,
    const std::vector<nova::fem::BoundaryEdge2D>& edges,
    const bool both_sides,
    const std::function<double(double, double)>& value
) {
    for (const auto& edge : edges) {
        const auto& a = problem.mesh.nodes[edge[0]];
        const auto& b = problem.mesh.nodes[edge[1]];
        if ((a.x == 0.0 && b.x == 0.0) ||
            (both_sides && a.x == 1.0 && b.x == 1.0)) {
            problem.dirichlet_edges.push_back({edge, value});
        }
    }
}

void require_exact_nodal_field(
    const nova::fem::ScalarDiffusion2DProblem& problem,
    nova::fem::H1TriangleOrder order,
    const nova::fem::ScalarDiffusion2DResult& result,
    const std::function<double(double, double)>& expected
) {
    NOVA_REQUIRE(result.linear_solver.converged);
    nova::fem::H1FunctionSpace2D space(problem.mesh, order);
    const auto exact_dofs = space.interpolate(expected);
    NOVA_REQUIRE(result.values.size() == space.dimension());
    for (std::size_t i = 0; i < exact_dofs.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.values[i], exact_dofs[i], 2.0e-9);
    }
}

} // namespace

NOVA_TEST(diffusion_h1_p2_reproduces_quadratic_manufactured_solution) {
    auto grid = nova::fem::make_rectangle_mesh_2d(0.0, 1.0, 0.0, 1.0, 3, 3);
    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.diffusivity = [](double, double) { return 1.0; };
    problem.source = [](double, double) { return 2.0; };
    const auto exact = [](double x, double) { return x * (1.0 - x); };
    constrain_vertical_boundaries(problem, grid.boundary_edges, true, exact);

    const auto result = nova::fem::solve_scalar_diffusion_h1_2d(
        problem, nova::fem::H1TriangleOrder::P2
    );
    NOVA_REQUIRE(result.values.size() > problem.mesh.nodes.size());
    require_exact_nodal_field(problem, nova::fem::H1TriangleOrder::P2, result, exact);
}

NOVA_TEST(diffusion_h1_p2_reproduces_linear_dirichlet_neumann_patch) {
    auto grid = nova::fem::make_rectangle_mesh_2d(0.0, 1.0, 0.0, 1.0, 3, 3);
    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.diffusivity = [](double, double) { return 1.0; };
    problem.source = [](double, double) { return 0.0; };
    constrain_vertical_boundaries(problem, grid.boundary_edges, false,
        [](double, double) { return 0.0; });

    for (const auto& edge : grid.boundary_edges) {
        const auto& a = problem.mesh.nodes[edge[0]];
        const auto& b = problem.mesh.nodes[edge[1]];
        if (a.x == 1.0 && b.x == 1.0) {
            problem.neumann_conditions.push_back({
                edge, [](double, double) { return -1.0; }
            });
        }
    }
    const auto result = nova::fem::solve_scalar_diffusion_h1_2d(
        problem, nova::fem::H1TriangleOrder::P2
    );
    require_exact_nodal_field(problem, nova::fem::H1TriangleOrder::P2, result,
        [](double x, double) { return x; });
}

NOVA_TEST(diffusion_h1_p2_reproduces_linear_dirichlet_robin_patch) {
    auto grid = nova::fem::make_rectangle_mesh_2d(0.0, 1.0, 0.0, 1.0, 3, 3);
    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.diffusivity = [](double, double) { return 1.0; };
    problem.source = [](double, double) { return 0.0; };
    constrain_vertical_boundaries(problem, grid.boundary_edges, false,
        [](double, double) { return 0.0; });

    for (const auto& edge : grid.boundary_edges) {
        const auto& a = problem.mesh.nodes[edge[0]];
        const auto& b = problem.mesh.nodes[edge[1]];
        if (a.x == 1.0 && b.x == 1.0) {
            problem.robin_conditions.push_back({
                edge,
                [](double, double) { return 2.0; },
                [](double, double) { return 1.5; }
            });
        }
    }
    const auto result = nova::fem::solve_scalar_diffusion_h1_2d(
        problem, nova::fem::H1TriangleOrder::P2
    );
    require_exact_nodal_field(problem, nova::fem::H1TriangleOrder::P2, result,
        [](double x, double) { return x; });
}

NOVA_TEST(diffusion_h1_p1_matches_existing_p1_on_linear_patch) {
    auto grid = nova::fem::make_rectangle_mesh_2d(0.0, 1.0, 0.0, 1.0, 4, 4);
    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.diffusivity = [](double, double) { return 1.0; };
    problem.source = [](double, double) { return 0.0; };
    for (std::size_t i : grid.boundary_nodes) {
        const auto& p = problem.mesh.nodes[i];
        problem.dirichlet_constraints.push_back({i, p.x + p.y});
    }
    const auto original = nova::fem::solve_scalar_diffusion_2d(problem);
    const auto generic = nova::fem::solve_scalar_diffusion_h1_2d(
        problem, nova::fem::H1TriangleOrder::P1
    );
    NOVA_REQUIRE(original.linear_solver.converged);
    NOVA_REQUIRE(generic.linear_solver.converged);
    NOVA_REQUIRE(original.values.size() == generic.values.size());
    for (std::size_t i = 0; i < generic.values.size(); ++i) {
        NOVA_REQUIRE_NEAR(generic.values[i], original.values[i], 1.0e-10);
    }
}

NOVA_TEST(diffusion_h1_p2_accepts_reversed_cell_orientation) {
    auto grid = nova::fem::make_rectangle_mesh_2d(0.0, 1.0, 0.0, 1.0, 2, 2);
    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    std::swap(problem.mesh.elements[0][1], problem.mesh.elements[0][2]);
    problem.diffusivity = [](double, double) { return 1.0; };
    problem.source = [](double, double) { return 2.0; };
    const auto exact = [](double x, double) { return x * (1.0 - x); };
    constrain_vertical_boundaries(problem, grid.boundary_edges, true, exact);
    const auto result = nova::fem::solve_scalar_diffusion_h1_2d(
        problem, nova::fem::H1TriangleOrder::P2
    );
    require_exact_nodal_field(problem, nova::fem::H1TriangleOrder::P2, result, exact);
}

NOVA_TEST(diffusion_h1_rejects_dirichlet_on_interior_edge) {
    auto grid = nova::fem::make_rectangle_mesh_2d(0.0, 1.0, 0.0, 1.0, 1, 1);
    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.diffusivity = [](double, double) { return 1.0; };
    problem.source = [](double, double) { return 0.0; };
    problem.dirichlet_edges.push_back({
        {0, 3}, [](double, double) { return 0.0; }
    });
    bool rejected = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_h1_2d(
            problem, nova::fem::H1TriangleOrder::P2
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}

NOVA_TEST(diffusion_h1_rejects_inconsistent_dirichlet_values) {
    auto grid = nova::fem::make_rectangle_mesh_2d(0.0, 1.0, 0.0, 1.0, 1, 1);
    nova::fem::ScalarDiffusion2DProblem problem;
    problem.mesh = std::move(grid.mesh);
    problem.diffusivity = [](double, double) { return 1.0; };
    problem.source = [](double, double) { return 0.0; };
    const auto edge = grid.boundary_edges.front();
    problem.dirichlet_constraints.push_back({edge[0], 1.0});
    problem.dirichlet_edges.push_back({
        edge, [](double, double) { return 0.0; }
    });
    bool rejected = false;
    try {
        static_cast<void>(nova::fem::solve_scalar_diffusion_h1_2d(
            problem, nova::fem::H1TriangleOrder::P2
        ));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}
