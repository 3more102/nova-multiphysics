#include "test_framework.hpp"
#include "nova/fem/diffusion_3d.hpp"

#include <cmath>
#include <vector>

namespace {
using namespace nova::fem;

bool face_plane(const TetrahedralMesh3D& mesh, const BoundaryFace3D& face,
                char axis, double coordinate) {
    for (const std::size_t vertex : face) {
        const auto& p = mesh.nodes[vertex];
        const double component = axis == 'x' ? p.x : (axis == 'y' ? p.y : p.z);
        if (component != coordinate) return false;
    }
    return true;
}

ScalarField3D constant(double value) {
    return [value](double, double, double) { return value; };
}

nova::linalg::SolverOptions accurate_options() {
    nova::linalg::SolverOptions options;
    options.preconditioner = nova::linalg::PreconditionerType::Jacobi;
    options.relative_tolerance = 1e-12;
    options.absolute_tolerance = 0.0;
    return options;
}
} // namespace

NOVA_TEST(energy_balance_3d_uses_uneliminated_dirichlet_reactions) {
    ScalarDiffusion3DProblem problem;
    problem.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    problem.diffusivity = constant(2.0);
    problem.source = constant(12.0);
    problem.cell_source = {6.0};
    problem.dirichlet_constraints = {{0,300}, {1,301}, {2,300}, {3,300}};

    // The default solve has no extra diagnostic storage or result.
    const auto ordinary = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(ordinary.linear_solver.converged);
    NOVA_REQUIRE(!ordinary.energy_balance.has_value());

    problem.compute_energy_balance = true;
    const auto result = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE(result.energy_balance.has_value());
    const auto& balance = *result.energy_balance;
    // Volume = 1/6 m3; source = 12 + 6 W/m3. With all nodes
    // prescribed, the reactions carry all 3 W of generated power.
    NOVA_REQUIRE_NEAR(balance.generated_heat_rate, 3.0, 2e-14);
    NOVA_REQUIRE_NEAR(balance.dirichlet_outward_heat_rate, 3.0, 2e-13);
    NOVA_REQUIRE_NEAR(balance.neumann_outward_heat_rate, 0.0, 0.0);
    NOVA_REQUIRE_NEAR(balance.robin_outward_heat_rate, 0.0, 0.0);
    NOVA_REQUIRE_NEAR(balance.net_outward_heat_rate, 3.0, 2e-13);
    NOVA_REQUIRE(std::abs(balance.imbalance) < 2e-13);
    NOVA_REQUIRE(balance.relative_imbalance < 1e-12);
}

NOVA_TEST(energy_balance_3d_closes_with_sources_flux_and_convection) {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 2, 2, 2);
    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = constant(2.0);
    problem.source = constant(20.0);
    problem.cell_source.assign(problem.mesh.elements.size(), 100.0);
    problem.compute_energy_balance = true;
    for (const auto& face : box.boundary_faces) {
        if (face_plane(problem.mesh, face, 'x', 0.0)) {
            problem.dirichlet_faces.push_back({face, constant(400)});
        } else if (face_plane(problem.mesh, face, 'x', 1.0)) {
            problem.robin_faces.push_back({face, constant(12), constant(300)});
        } else if (face_plane(problem.mesh, face, 'y', 1.0)) {
            problem.neumann_faces.push_back({face, constant(20)});
        }
    }

    const auto result = solve_scalar_diffusion_3d(problem, accurate_options());
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE(result.energy_balance.has_value());
    const auto& balance = *result.energy_balance;
    NOVA_REQUIRE_NEAR(balance.generated_heat_rate, 120.0, 1e-10);
    NOVA_REQUIRE_NEAR(balance.neumann_outward_heat_rate, 20.0, 1e-11);
    NOVA_REQUIRE(balance.robin_outward_heat_rate > 0.0);
    NOVA_REQUIRE(balance.dirichlet_outward_heat_rate < 0.0);
    NOVA_REQUIRE_NEAR(balance.net_outward_heat_rate, 120.0, 2e-8);
    NOVA_REQUIRE(std::abs(balance.imbalance) < 2e-8);
    NOVA_REQUIRE(balance.relative_imbalance < 1e-9);
}

NOVA_TEST(energy_balance_3d_neumann_nullspace_gauge_is_not_physical_reaction) {
    const auto box = make_box_mesh_3d(0, 1, 0, 1, 0, 1, 1, 1, 1);
    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = constant(2.0);
    problem.source = constant(0.0);
    problem.compute_energy_balance = true;
    for (const auto& face : box.boundary_faces) {
        if (face_plane(problem.mesh, face, 'x', 0.0))
            problem.neumann_faces.push_back({face, constant(-5.0)});
        if (face_plane(problem.mesh, face, 'x', 1.0))
            problem.neumann_faces.push_back({face, constant(5.0)});
    }
    const auto result = solve_scalar_diffusion_3d(problem, accurate_options());
    NOVA_REQUIRE(result.linear_solver.converged);
    NOVA_REQUIRE(result.energy_balance.has_value());
    const auto& balance = *result.energy_balance;
    NOVA_REQUIRE_NEAR(balance.dirichlet_outward_heat_rate, 0.0, 0.0);
    NOVA_REQUIRE_NEAR(balance.generated_heat_rate, 0.0, 0.0);
    NOVA_REQUIRE_NEAR(balance.neumann_outward_heat_rate, 0.0, 1e-13);
    NOVA_REQUIRE_NEAR(balance.imbalance, 0.0, 1e-13);
}
