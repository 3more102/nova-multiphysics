#include "test_framework.hpp"

#include "nova/fem/diffusion_2d.hpp"
#include "nova/fem/h1_space_2d.hpp"

#include <cmath>
#include <limits>
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

NOVA_TEST(neumann_2d_reproduces_affine_fields_with_geometric_zero_mean) {
    using namespace nova::fem;
    auto grid = make_rectangle_mesh_2d(0, 1, 0, 1, 2, 2);
    ScalarDiffusion2DProblem p;
    p.mesh = grid.mesh;
    for (auto& n : p.mesh.nodes) {
        if (n.x == 0.5) n.x = 0.13;
    }
    std::swap(p.mesh.elements[1][0], p.mesh.elements[1][1]);
    p.diffusivity = [](double, double) { return 1.0; };
    p.source = [](double, double) { return 0.0; };
    for (const auto& e : grid.boundary_edges) {
        const auto& a = p.mesh.nodes[e[0]];
        const auto& b = p.mesh.nodes[e[1]];
        if (a.x == b.x) {
            const double q = a.x == 0.0 ? 1.0 : -1.0;
            p.neumann_conditions.push_back({{e[1], e[0]}, [q](double, double) { return q; }});
        } else {
            p.robin_conditions.push_back({e, p.source, p.diffusivity}); // zero h is insulating
        }
    }
    for (auto order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        const auto r = solve_scalar_diffusion_h1_2d(p, order);
        require_exact_nodal_field(p, order, r, [](double x, double) { return x - 0.5; });
        NOVA_REQUIRE_NEAR(r.linear_solver.solution.front(), 0.0, 1e-13);
        if (order == H1TriangleOrder::P1) {
            NOVA_REQUIRE(solve_scalar_diffusion_2d(p).values == r.values);
        }
    }
}

NOVA_TEST(neumann_p2_uses_integrated_quadratic_field_and_balanced_volume_source) {
    using namespace nova::fem;
    auto grid = make_rectangle_mesh_2d(0, 1, 0, 1, 3, 2);
    ScalarDiffusion2DProblem p;
    p.mesh = grid.mesh;
    p.diffusivity = [](double, double) { return 1.0; };
    p.source = [](double, double) { return -2.0; };
    for (const auto& e : grid.boundary_edges) {
        if (p.mesh.nodes[e[0]].x == 1.0 && p.mesh.nodes[e[1]].x == 1.0) {
            p.neumann_conditions.push_back({e, p.source});
        }
    }
    const auto r = solve_scalar_diffusion_h1_2d(p, H1TriangleOrder::P2);
    require_exact_nodal_field(p, H1TriangleOrder::P2, r,
        [](double x, double) { return x*x - 1.0/3.0; });
}

NOVA_TEST(neumann_2d_handles_floating_robin_and_dirichlet_components_in_one_solve) {
    using namespace nova::fem;
    ScalarDiffusion2DProblem p;
    p.mesh = {{{0,0},{1,0},{0,1}, {2,0},{3,0},{2,1}, {4,0},{5,0},{4,1}},
              {{0,1,2}, {3,5,4}, {6,7,8}}};
    p.diffusivity = [](double, double) { return 1.0; };
    p.source = [](double, double) { return 0.0; };
    p.neumann_conditions = {{{0,2}, p.diffusivity},
        {{1,2}, [](double, double) { return -1.0/std::sqrt(2.0); }}};
    p.robin_conditions = {{{4,5}, p.diffusivity, [](double, double) { return 5.0; }}};
    p.dirichlet_constraints.push_back({6,9.0});
    for (auto order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        const auto r = solve_scalar_diffusion_h1_2d(p, order);
        require_exact_nodal_field(p, order, r, [](double x, double) {
            if (x < 2.0) return x - 1.0/3.0;
            return x < 4.0 ? 5.0 : 9.0;
        });
    }
}

NOVA_TEST(neumann_2d_rejects_global_cancellation_and_tiny_unbalanced_components) {
    using namespace nova::fem;
    ScalarDiffusion2DProblem p;
    p.mesh = {{{0,0},{1,0},{0,1}, {2,0},{3,0},{2,1}}, {{0,1,2}, {3,4,5}}};
    p.diffusivity = [](double, double) { return 1.0; };
    for (double magnitude : {1.0, 1.0e-24}) {
        p.source = [magnitude](double x, double) { return x < 2 ? magnitude : -magnitude; };
        for (auto order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
            for (bool anchor_first : {false, true}) {
                p.dirichlet_constraints.clear();
                if (anchor_first) p.dirichlet_constraints.push_back({0,1.0});
                bool rejected = false;
                try {
                    (void)solve_scalar_diffusion_h1_2d(p, order);
                } catch (const std::invalid_argument&) {
                    rejected = true;
                }
                NOVA_REQUIRE(rejected);
            }
        }
    }
}

NOVA_TEST(stationary_diffusion_preserves_prescribed_identity_rows_at_solver_tolerance) {
    using namespace nova::fem;
    auto grid = make_rectangle_mesh_2d(0, 1, 0, 1, 4, 4);
    ScalarDiffusion2DProblem p;
    p.mesh = grid.mesh;
    p.diffusivity = [](double, double) { return 1.0; };
    p.source = [](double x, double y) { return x + y; };
    p.dirichlet_constraints = {{0,3.25}, {17,-1.75}};
    nova::linalg::SolverOptions options;
    options.relative_tolerance = 0.05;
    for (auto order : {H1TriangleOrder::P1, H1TriangleOrder::P2}) {
        const auto r = solve_scalar_diffusion_h1_2d(p, order, options);
        NOVA_REQUIRE(r.linear_solver.converged);
        for (const auto& boundary : p.dirichlet_constraints) {
            NOVA_REQUIRE(r.values[boundary.node] == boundary.value);
            NOVA_REQUIRE(r.linear_solver.solution[boundary.node] == boundary.value);
        }
    }
}

NOVA_TEST(diffusion_rejects_overflow_even_when_all_stiffness_rows_are_prescribed) {
    using namespace nova::fem;
    ScalarDiffusion2DProblem p;
    p.mesh = {{{0,0}, {1,0}, {0,1.0e-10}}, {{0,1,2}}};
    p.diffusivity = [](double, double) { return std::numeric_limits<double>::max(); };
    p.source = [](double, double) { return 0.0; };
    p.dirichlet_constraints = {{0,0.0}, {1,0.0}, {2,0.0}};
    bool rejected = false;
    try {
        (void)solve_scalar_diffusion_2d(p);
    } catch (const std::overflow_error&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}
