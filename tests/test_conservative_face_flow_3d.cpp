#include "test_framework.hpp"

#include "nova/fem/conservative_face_flow_3d.hpp"
#include "nova/fem/energy_balance_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using namespace nova::fem;

ScalarDiffusion3DProblem one_tetra(double source) {
    ScalarDiffusion3DProblem problem;
    problem.mesh = {{{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}, {{0,1,2,3}}};
    problem.diffusivity = [](double, double, double) { return 2.0; };
    problem.source = [source](double, double, double) { return source; };
    return problem;
}

void check_cellwise_conservation(const ConservativeHeatFlow3D& result) {
    std::vector<double> net(result.generated_heat_per_cell.size(), 0.0);
    for (const auto& face : result.faces) {
        net[face.owner_cell] += face.outward_heat_rate;
        if (face.neighbor_cell) net[*face.neighbor_cell] -= face.outward_heat_rate;
    }
    for (std::size_t cell = 0; cell < net.size(); ++cell) {
        NOVA_REQUIRE_NEAR(net[cell], result.generated_heat_per_cell[cell], 2e-10);
    }
    NOVA_REQUIRE(result.maximum_defect_after <= 2e-10);
}

template<class F>
bool rejects(F&& f) {
    try { f(); }
    catch (const std::invalid_argument&) { return true; }
    return false;
}

} // namespace

NOVA_TEST(conservative_face_flow_assigns_generated_heat_to_dirichlet_outlet) {
    auto problem = one_tetra(6.0); // Volume 1/6; 1 W total.
    problem.dirichlet_faces.push_back({
        {0,1,2}, [](double, double, double) { return 0.0; }
    });
    const auto result = reconstruct_conservative_heat_flow_3d(
        problem, std::array{0.0,0.0,0.0,0.0});
    NOVA_REQUIRE(result.faces.size() == 4);
    NOVA_REQUIRE_NEAR(result.total_generated_heat, 1.0, 1e-14);
    NOVA_REQUIRE_NEAR(result.total_exterior_outward_heat, 1.0, 1e-14);
    NOVA_REQUIRE_NEAR(result.maximum_defect_before, 1.0, 1e-14);
    NOVA_REQUIRE_NEAR(result.maximum_absolute_face_correction, 1.0, 1e-14);
    std::size_t outlets = 0;
    for (const auto& face : result.faces) {
        if (face.kind == FaceFlowBoundaryKind3D::dirichlet) {
            ++outlets;
            NOVA_REQUIRE_NEAR(face.outward_heat_rate, 1.0, 1e-14);
        } else {
            NOVA_REQUIRE(face.kind == FaceFlowBoundaryKind3D::insulated);
            NOVA_REQUIRE_NEAR(face.outward_heat_rate, 0.0, 1e-14);
        }
    }
    NOVA_REQUIRE(outlets == 1);
    check_cellwise_conservation(result);
}

NOVA_TEST(conservative_face_flow_preserves_exact_neumann_heat_without_outlet) {
    auto problem = one_tetra(6.0);
    problem.neumann_faces.push_back({
        {0,2,3}, [](double, double, double) { return 2.0; }
    }); // Area = 1/2, so 1 W leaves.
    const auto result = reconstruct_conservative_heat_flow_3d(
        problem, std::array{0.0,0.0,0.0,0.0});
    NOVA_REQUIRE_NEAR(result.total_exterior_outward_heat, 1.0, 1e-14);
    NOVA_REQUIRE_NEAR(result.maximum_absolute_face_correction, 0.0, 1e-14);
    for (const auto& face : result.faces) {
        if (face.kind == FaceFlowBoundaryKind3D::neumann)
            NOVA_REQUIRE_NEAR(face.outward_heat_rate, 1.0, 1e-14);
        if (face.kind == FaceFlowBoundaryKind3D::insulated)
            NOVA_REQUIRE_NEAR(face.outward_heat_rate, 0.0, 1e-14);
    }
    check_cellwise_conservation(result);
}

NOVA_TEST(conservative_face_flow_has_one_shared_rate_for_discontinuous_anisotropic_materials) {
    auto problem = one_tetra(0.0);
    problem.mesh.nodes.push_back({1,1,1});
    problem.mesh.elements.push_back({4,1,3,2});
    problem.diffusivity = {};
    problem.cell_diffusivity_tensors = {
        {1.0,2.0,3.0,0.0,0.0,0.0},
        {4.0,5.0,6.0,0.0,0.0,0.0}
    };
    problem.cell_source = {6.0, 12.0};
    problem.dirichlet_faces.push_back({
        {0,1,2}, [](double, double, double) { return 0.0; }
    });
    const auto result = reconstruct_conservative_heat_flow_3d(
        problem, std::array{0.0,0.0,0.0,0.0,0.0});
    NOVA_REQUIRE(result.faces.size() == 7);
    std::size_t shared = 0;
    for (const auto& face : result.faces) {
        if (face.neighbor_cell) {
            ++shared;
            NOVA_REQUIRE(face.owner_cell == 0);
            NOVA_REQUIRE(*face.neighbor_cell == 1);
            NOVA_REQUIRE(std::abs(face.outward_heat_rate) > 0.0);
        }
    }
    NOVA_REQUIRE(shared == 1);
    NOVA_REQUIRE_NEAR(result.total_exterior_outward_heat,
                      result.total_generated_heat, 1e-12);
    check_cellwise_conservation(result);
}

NOVA_TEST(conservative_face_flow_respects_robin_faces_and_tagged_sources) {
    auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    ScalarDiffusion3DProblem problem;
    problem.mesh = std::move(box.mesh);
    problem.source = [](double, double, double) { return 0.0; };
    for (const auto& cell : problem.mesh.elements) {
        double mean_x = 0.0;
        for (auto node : cell) mean_x += problem.mesh.nodes[node].x / 4.0;
        const double k = mean_x < 0.5 ? 2.0 : 6.0;
        problem.cell_diffusivity_tensors.push_back({k, 4.0, 5.0, 0.0, 0.0, 0.0});
        problem.cell_source.push_back(mean_x < 0.5 ? 120.0 : 360.0);
    }
    for (const auto& face : box.boundary_faces) {
        const auto on_x = [&](double x) {
            return std::all_of(face.begin(), face.end(),
                [&](std::size_t node) { return problem.mesh.nodes[node].x == x; });
        };
        if (on_x(0.0)) {
            problem.dirichlet_faces.push_back({
                face, [](double, double, double) { return 400.0; }
            });
        } else if (on_x(1.0)) {
            problem.robin_faces.push_back({
                face, [](double, double, double) { return 12.0; },
                [](double, double, double) { return 300.0; }
            });
        }
    }
    const auto solved = solve_scalar_diffusion_3d(problem);
    NOVA_REQUIRE(solved.linear_solver.converged);
    const auto weak = evaluate_diffusion_energy_balance_3d(problem, solved.values);
    NOVA_REQUIRE(std::abs(weak.signed_balance_defect) < 2e-8);
    const auto result = reconstruct_conservative_heat_flow_3d(problem, solved.values);
    NOVA_REQUIRE_NEAR(result.total_generated_heat, 240.0, 2e-9);
    NOVA_REQUIRE_NEAR(result.total_exterior_outward_heat, 240.0, 2e-8);
    long double robin = 0.0L;
    for (const auto& face : result.faces) {
        if (face.kind == FaceFlowBoundaryKind3D::robin) {
            robin += face.outward_heat_rate;
            NOVA_REQUIRE_NEAR(face.outward_heat_rate,
                              face.raw_outward_heat_rate, 1e-12);
        }
        if (face.kind == FaceFlowBoundaryKind3D::insulated) {
            NOVA_REQUIRE_NEAR(face.outward_heat_rate, 0.0, 1e-14);
        }
    }
    NOVA_REQUIRE_NEAR(robin, weak.robin_outward_heat, 2e-8);
    check_cellwise_conservation(result);
}


NOVA_TEST(conservative_face_flow_balances_each_disconnected_component_separately) {
    auto problem = one_tetra(6.0);
    problem.mesh.nodes.insert(problem.mesh.nodes.end(), {
        {3,0,0}, {4,0,0}, {3,1,0}, {3,0,1}
    });
    problem.mesh.elements.push_back({4,5,6,7});
    problem.dirichlet_faces.push_back({
        {0,1,2}, [](double, double, double) { return 0.0; }
    });
    problem.dirichlet_faces.push_back({
        {4,5,6}, [](double, double, double) { return 0.0; }
    });
    const std::array zeros{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    const auto balanced = reconstruct_conservative_heat_flow_3d(problem, zeros);
    NOVA_REQUIRE_NEAR(balanced.total_generated_heat, 2.0, 1e-14);
    NOVA_REQUIRE_NEAR(balanced.total_exterior_outward_heat, 2.0, 1e-14);
    check_cellwise_conservation(balanced);
    // The other component cannot discharge its source through an outlet
    // owned by the first disconnected component.
    problem.dirichlet_faces.pop_back();
    NOVA_REQUIRE(rejects([&] {
        (void)reconstruct_conservative_heat_flow_3d(problem, zeros);
    }));
}

NOVA_TEST(conservative_face_flow_rejects_incompatible_or_invalid_inputs) {
    auto problem = one_tetra(6.0);
    const std::array zero{0.0,0.0,0.0,0.0};
    NOVA_REQUIRE(rejects([&] {
        (void)reconstruct_conservative_heat_flow_3d(problem, zero);
    }));
    problem.dirichlet_faces.push_back({
        {0,1,2}, [](double, double, double) { return 0.0; }
    });
    NOVA_REQUIRE(rejects([&] {
        (void)reconstruct_conservative_heat_flow_3d(problem,
                                                    std::array{0.0,1.0,0.0,0.0});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)reconstruct_conservative_heat_flow_3d(
            problem, std::array{0.0,0.0,std::numeric_limits<double>::quiet_NaN(),0.0});
    }));
    NOVA_REQUIRE(rejects([&] {
        (void)reconstruct_conservative_heat_flow_3d(problem, std::array{0.0,0.0});
    }));
    problem.neumann_faces.push_back({
        {0,1,2}, [](double, double, double) { return 1.0; }
    });
    NOVA_REQUIRE(rejects([&] {
        (void)reconstruct_conservative_heat_flow_3d(problem, zero);
    }));
}
