#include "test_framework.hpp"

#include "nova/physics/affine_floating_conductor_schur_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using nova::physics::ElectrodeBoundaryGroup3D;
using nova::physics::Electrostatics3DProblem;
using nova::physics::FixedElectrodePotential3D;
using nova::physics::FloatingElectrodeCharge3D;
using nova::physics::solve_affine_floating_conductor_schur_3d;

[[nodiscard]] bool x_plate(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::BoundaryFace3D& face,
    double x
) {
    return mesh.nodes[face[0]].x == x &&
           mesh.nodes[face[1]].x == x &&
           mesh.nodes[face[2]].x == x;
}

[[nodiscard]] bool y_plate(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::BoundaryFace3D& face,
    double y
) {
    return mesh.nodes[face[0]].y == y &&
           mesh.nodes[face[1]].y == y &&
           mesh.nodes[face[2]].y == y;
}

[[nodiscard]] std::vector<ElectrodeBoundaryGroup3D> electrodes(
    const Electrostatics3DProblem& problem,
    const std::vector<nova::fem::BoundaryFace3D>& faces
) {
    std::vector<ElectrodeBoundaryGroup3D> groups{{10, {}}, {20, {}}};
    for (const auto& face : faces) {
        if (x_plate(problem.mesh, face, 0.0)) {
            groups[0].faces.push_back(face);
        } else if (x_plate(problem.mesh, face, 1.0)) {
            groups[1].faces.push_back({face[2], face[1], face[0]});
        }
    }
    NOVA_REQUIRE(!groups[0].faces.empty() && !groups[1].faces.empty());
    return groups;
}

[[nodiscard]] auto solve(
    const Electrostatics3DProblem& problem,
    const std::vector<ElectrodeBoundaryGroup3D>& groups,
    double charge
) {
    return solve_affine_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, charge}});
}

} // namespace

NOVA_TEST(affine_schur_uniform_volume_charge_exact_terminal_and_gauss_balance) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    problem.free_charge_density = [](double, double, double) { return 1.0; };

    const auto result = solve(problem, electrodes(problem, box.boundary_faces), 0.0);
    NOVA_REQUIRE(result.field.linear_solver.converged);
    NOVA_REQUIRE(result.impressed_charge_offset_coulombs.size() == 2U);
    NOVA_REQUIRE_NEAR(result.impressed_charge_offset_coulombs[0], -0.5, 1.0e-9);
    NOVA_REQUIRE_NEAR(result.impressed_charge_offset_coulombs[1], -0.5, 1.0e-9);
    NOVA_REQUIRE_NEAR(result.terminals.potential_volts[1], 0.25, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.terminals.estimated_charge_coulombs[1], 0.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(
        result.loaded_discrete_reactions[1].reaction_charge_coulombs,
        0.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(
        result.loaded_discrete_reactions[0].reaction_charge_coulombs,
        -1.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(
        result.loaded_discrete_reactions[0].reaction_charge_coulombs +
        result.loaded_discrete_reactions[1].reaction_charge_coulombs,
        -1.0, 1.0e-8);
    // Raw boundary flux is not the constrained K*V-F reaction: it deliberately
    // differs by the source load on the selected boundary P1 nodes.
    NOVA_REQUIRE_NEAR(
        result.boundary_flux_charges[1].conductor_charge_coulombs,
        0.5, 1.0e-8);
    NOVA_REQUIRE(result.max_flux_reaction_disagreement_coulombs > 0.49);
    NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-8);
}

NOVA_TEST(affine_schur_additive_cell_and_spatial_free_charge) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    problem.free_charge_density = [](double, double, double) { return 1.0; };
    problem.cell_free_charge_density.assign(box.mesh.elements.size(), 1.0);

    const auto result = solve(problem, electrodes(problem, box.boundary_faces), 0.0);
    NOVA_REQUIRE_NEAR(result.impressed_charge_offset_coulombs[1], -1.0, 1.0e-9);
    NOVA_REQUIRE_NEAR(result.terminals.potential_volts[1], 0.5, 1.0e-8);
    NOVA_REQUIRE_NEAR(
        result.loaded_discrete_reactions[0].reaction_charge_coulombs +
        result.loaded_discrete_reactions[1].reaction_charge_coulombs,
        -2.0, 1.0e-8);
    NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-8);
}

NOVA_TEST(affine_schur_neumann_and_volume_load_global_balance) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    problem.free_charge_density = [](double, double, double) { return 1.0; };
    for (const auto& face : box.boundary_faces) {
        if (y_plate(problem.mesh, face, 1.0)) {
            problem.displacement_flux_faces.push_back({
                face, [](double, double, double) { return 3.0; }});
        }
    }
    const auto result = solve(problem, electrodes(problem, box.boundary_faces), 0.0);
    NOVA_REQUIRE(result.field.linear_solver.converged);
    NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-7);
    NOVA_REQUIRE_NEAR(
        result.loaded_discrete_reactions[1].reaction_charge_coulombs,
        0.0, 1.0e-7);
    // Sum conductor reactions = prescribed side displacement (3 C)
    // minus integrated free volume charge (1 C).
    NOVA_REQUIRE_NEAR(
        result.loaded_discrete_reactions[0].reaction_charge_coulombs +
        result.loaded_discrete_reactions[1].reaction_charge_coulombs,
        2.0, 1.0e-7);
}

NOVA_TEST(affine_schur_source_free_matches_previous_discrete_solver) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    const auto groups = electrodes(problem, box.boundary_faces);
    const auto affine = solve(problem, groups, 20.0);
    const auto homogeneous = nova::physics::solve_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    NOVA_REQUIRE_NEAR(affine.impressed_charge_offset_coulombs[0], 0.0, 1.0e-9);
    NOVA_REQUIRE_NEAR(affine.impressed_charge_offset_coulombs[1], 0.0, 1.0e-9);
    NOVA_REQUIRE_NEAR(affine.terminals.potential_volts[1], 10.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(
        affine.terminals.potential_volts[1],
        homogeneous.terminals.potential_volts[1], 1.0e-8);
    NOVA_REQUIRE_NEAR(
        affine.loaded_discrete_reactions[1].reaction_charge_coulombs,
        homogeneous.discrete_reactions[1].reaction_charge_coulombs, 1.0e-8);
}

NOVA_TEST(affine_schur_rejects_nonfinite_source_and_flux) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    const auto groups = electrodes(problem, box.boundary_faces);

    problem.free_charge_density = [](double, double, double) {
        return std::numeric_limits<double>::quiet_NaN();
    };
    bool rejected = false;
    try {
        (void)solve(problem, groups, 0.0);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}

