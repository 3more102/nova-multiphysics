#include "test_framework.hpp"

#include "nova/physics/floating_conductor_schur_3d.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

using nova::physics::Electrostatics3DProblem;
using nova::physics::ElectrodeBoundaryGroup3D;
using nova::physics::FixedElectrodePotential3D;
using nova::physics::FloatingElectrodeCharge3D;
using nova::physics::solve_floating_conductor_schur_3d;

[[nodiscard]] bool on_x(
    const nova::fem::TetrahedralMesh3D& mesh,
    const nova::fem::BoundaryFace3D& face,
    const double x
) {
    return mesh.nodes[face[0]].x == x &&
           mesh.nodes[face[1]].x == x &&
           mesh.nodes[face[2]].x == x;
}

[[nodiscard]] std::vector<ElectrodeBoundaryGroup3D> plates(
    const Electrostatics3DProblem& problem,
    const std::vector<nova::fem::BoundaryFace3D>& boundary_faces
) {
    std::vector<ElectrodeBoundaryGroup3D> groups{{10, {}}, {20, {}}};
    for (const auto& face : boundary_faces) {
        if (on_x(problem.mesh, face, 0.0)) {
            groups[0].faces.push_back(face);
        } else if (on_x(problem.mesh, face, 1.0)) {
            groups[1].faces.push_back({face[2], face[1], face[0]});
        }
    }
    NOVA_REQUIRE(!groups[0].faces.empty() && !groups[1].faces.empty());
    return groups;
}

template <class Callable>
void expect_invalid(const Callable& operation) {
    bool invalid = false;
    try {
        operation();
    } catch (const std::invalid_argument&) {
        invalid = true;
    }
    NOVA_REQUIRE(invalid);
}

} // namespace

NOVA_TEST(floating_schur_parallel_plate_exact_charge_voltage_energy) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    const auto groups = plates(problem, box.boundary_faces);

    const auto result = solve_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    NOVA_REQUIRE(result.field.linear_solver.converged);
    NOVA_REQUIRE(result.discrete_capacitance.electrode_ids ==
                 std::vector<int>({10, 20}));
    NOVA_REQUIRE(result.discrete_capacitance.matrix_farads.size() == 4U);
    NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[0], 2.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[1], -2.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[2], -2.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[3], 2.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.terminals.potential_volts[0], 0.0, 1.0e-9);
    NOVA_REQUIRE_NEAR(result.terminals.potential_volts[1], 10.0, 1.0e-9);
    NOVA_REQUIRE(result.discrete_reactions.size() == 2U);
    NOVA_REQUIRE(result.boundary_flux_charges.size() == 2U);
    NOVA_REQUIRE_NEAR(result.discrete_reactions[0].reaction_charge_coulombs,
                      -20.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.discrete_reactions[1].reaction_charge_coulombs,
                      20.0, 1.0e-8);
    NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-8);
    NOVA_REQUIRE(result.max_flux_reaction_disagreement_coulombs < 1.0e-8);
    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.field.potential[i],
                          10.0 * problem.mesh.nodes[i].x, 1.0e-8);
    }
    for (const auto& field : result.field.electric_field) {
        NOVA_REQUIRE_NEAR(field.x, -10.0, 1.0e-8);
        NOVA_REQUIRE_NEAR(field.y, 0.0, 1.0e-8);
        NOVA_REQUIRE_NEAR(field.z, 0.0, 1.0e-8);
    }
    NOVA_REQUIRE_NEAR(result.field.stored_energy_joules, 100.0, 1.0e-7);
    NOVA_REQUIRE_NEAR(
        result.terminals.potential_volts[1] *
            result.discrete_reactions[1].reaction_charge_coulombs,
        2.0 * result.field.stored_energy_joules, 1.0e-7);
}

NOVA_TEST(floating_schur_variable_permittivity_charge_and_energy_identity) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 3U, 2U, 2U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double x, double y, double) {
        return 1.0 + x*x + 0.25*y;
    };
    const auto groups = plates(problem, box.boundary_faces);
    const auto result = solve_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, 7.0}});
    NOVA_REQUIRE(result.field.linear_solver.converged);
    NOVA_REQUIRE(result.terminals.potential_volts[1] > 0.0);
    NOVA_REQUIRE(result.discrete_reactions.size() == 2U);
    NOVA_REQUIRE_NEAR(result.discrete_reactions[1].reaction_charge_coulombs,
                      7.0, 2.0e-7);
    NOVA_REQUIRE_NEAR(result.discrete_reactions[0].reaction_charge_coulombs,
                      -7.0, 2.0e-7);
    NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 2.0e-7);
    NOVA_REQUIRE(result.discrete_capacitance.max_reciprocity_defect_farads <
                 2.0e-8);
    NOVA_REQUIRE_NEAR(
        result.terminals.potential_volts[1] *
            result.discrete_reactions[1].reaction_charge_coulombs,
        2.0 * result.field.stored_energy_joules, 2.0e-6);
}

NOVA_TEST(floating_schur_single_floating_with_external_ground) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    ElectrodeBoundaryGroup3D high{20, {}};
    for (const auto& face : box.boundary_faces) {
        if (on_x(problem.mesh, face, 0.0)) {
            // An existing non-electrode potential is deliberately nonzero.
            // The capacitance extraction interprets it as a 0 V reference.
            problem.potential_faces.push_back({
                face, [](double, double, double) { return 9.0; }});
        } else if (on_x(problem.mesh, face, 1.0)) {
            high.faces.push_back(face);
        }
    }
    const auto result = solve_floating_conductor_schur_3d(
        problem, std::vector<ElectrodeBoundaryGroup3D>{high},
        std::vector<FixedElectrodePotential3D>{},
        std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    NOVA_REQUIRE(result.field.linear_solver.converged);
    NOVA_REQUIRE(result.discrete_capacitance.matrix_farads.size() == 1U);
    NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[0], 2.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.terminals.potential_volts[0], 10.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.discrete_reactions[0].reaction_charge_coulombs,
                      20.0, 1.0e-8);
    NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-8);
    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.field.potential[i],
                          10.0 * problem.mesh.nodes[i].x, 1.0e-8);
    }
}

NOVA_TEST(floating_schur_all_fixed_terminal_voltages) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    const auto result = solve_floating_conductor_schur_3d(
        problem, plates(problem, box.boundary_faces),
        std::vector<FixedElectrodePotential3D>{{20, 10.0}, {10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{});
    NOVA_REQUIRE(result.field.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.terminals.potential_volts[1], 10.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.discrete_reactions[0].reaction_charge_coulombs,
                      -20.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.discrete_reactions[1].reaction_charge_coulombs,
                      20.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(result.max_floating_reaction_mismatch_coulombs, 0.0,
                      1.0e-12);
}

NOVA_TEST(floating_schur_rejects_invalid_loads_gauge_and_ids) {
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    const auto groups = plates(problem, box.boundary_faces);
    const auto check = [&](const std::vector<FixedElectrodePotential3D>& fixed,
                           const std::vector<FloatingElectrodeCharge3D>& floating) {
        (void)solve_floating_conductor_schur_3d(problem, groups, fixed, floating);
    };
    // Both electrodes left unreferenced: Maxwell capacitance has a gauge mode.
    expect_invalid([&] {
        check({}, {{10, -20.0}, {20, 20.0}});
    });
    expect_invalid([&] { check({{10, 0.0}}, {{10, 20.0}}); });
    expect_invalid([&] { check({{99, 0.0}}, {{20, 20.0}}); });
    expect_invalid([&] { check({{10, 0.0}}, {}); });
    expect_invalid([&] {
        check({{10, 0.0}, {10, 0.0}}, {{20, 20.0}});
    });

    auto loaded = problem;
    loaded.free_charge_density = [](double, double, double) { return 1.0; };
    expect_invalid([&] {
        (void)solve_floating_conductor_schur_3d(
            loaded, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
            std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    });
    loaded.free_charge_density = {};
    loaded.cell_free_charge_density.assign(problem.mesh.elements.size(), 1.0);
    expect_invalid([&] {
        (void)solve_floating_conductor_schur_3d(
            loaded, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
            std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    });
    loaded.cell_free_charge_density.clear();
    loaded.displacement_flux_faces.push_back({
        box.boundary_faces.front(),
        [](double, double, double) { return 1.0; }
    });
    expect_invalid([&] {
        (void)solve_floating_conductor_schur_3d(
            loaded, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
            std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    });
}


NOVA_TEST(floating_schur_anisotropic_tensor_gauge_shift_patch) {
    // On a unit cube, a diagonal SPD permittivity tensor with
    // epsilon_xx=2 F/m and insulated y/z sidewalls gives the exact field
    // V(x)=V_left+3x when the right plate carries +6 C.
    // Both the callback and the cellwise tensor assembly must return
    // Q_left=-6 C, Q_right=+6 C, E_x=-3 V/m, D_x=-6 C/m^2,
    // U=9 J, and the same C=[[2,-2],[-2,2]] F.
    //
    // Shifting the voltage reference changes neither E, D, charges,
    // capacitance nor stored energy, while every nodal V shifts equally.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    for (int material = 0; material != 2; ++material) {
        Electrostatics3DProblem problem;
        problem.mesh = box.mesh;
        if (material == 0) {
            problem.permittivity_tensor = [](double, double, double) {
                return nova::fem::SymmetricTensor3D{2.0, 5.0, 7.0, 0.0, 0.0, 0.0};
            };
        } else {
            problem.cell_permittivity_tensors.assign(
                problem.mesh.elements.size(),
                nova::fem::SymmetricTensor3D{2.0, 5.0, 7.0, 0.0, 0.0, 0.0});
        }
        const auto groups = plates(problem, box.boundary_faces);
        for (const double reference : std::vector<double>{-3.0, 4.0}) {
            const auto result = solve_floating_conductor_schur_3d(
                problem, groups,
                std::vector<FixedElectrodePotential3D>{{10, reference}},
                std::vector<FloatingElectrodeCharge3D>{{20, 6.0}});
            NOVA_REQUIRE(result.field.linear_solver.converged);
            NOVA_REQUIRE(result.discrete_capacitance.matrix_farads.size() == 4U);
            NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[0],
                              2.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[1],
                              -2.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[2],
                              -2.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(result.discrete_capacitance.matrix_farads[3],
                              2.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(result.terminals.potential_volts[0],
                              reference, 1.0e-8);
            NOVA_REQUIRE_NEAR(result.terminals.potential_volts[1],
                              reference + 3.0, 1.0e-8);
            NOVA_REQUIRE(result.discrete_reactions.size() == 2U);
            NOVA_REQUIRE_NEAR(result.discrete_reactions[0].reaction_charge_coulombs,
                              -6.0, 1.0e-8);
            NOVA_REQUIRE_NEAR(result.discrete_reactions[1].reaction_charge_coulombs,
                              6.0, 1.0e-8);
            NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-8);
            NOVA_REQUIRE(result.max_flux_reaction_disagreement_coulombs < 1.0e-8);
            NOVA_REQUIRE_NEAR(result.field.stored_energy_joules, 9.0, 1.0e-8);
            for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
                NOVA_REQUIRE_NEAR(result.field.potential[i],
                    reference + 3.0 * problem.mesh.nodes[i].x, 1.0e-8);
            }
            for (std::size_t i = 0; i < problem.mesh.elements.size(); ++i) {
                const auto& e = result.field.electric_field[i];
                const auto& d = result.field.displacement[i];
                NOVA_REQUIRE_NEAR(e.x, -3.0, 1.0e-8);
                NOVA_REQUIRE_NEAR(e.y, 0.0, 1.0e-8);
                NOVA_REQUIRE_NEAR(e.z, 0.0, 1.0e-8);
                NOVA_REQUIRE_NEAR(d.x, -6.0, 1.0e-8);
                NOVA_REQUIRE_NEAR(d.y, 0.0, 1.0e-8);
                NOVA_REQUIRE_NEAR(d.z, 0.0, 1.0e-8);
            }
            // Gauge-invariant electrostatic energy, with nonzero fixed voltage.
            const double terminal_work =
                result.terminals.potential_volts[0] *
                    result.discrete_reactions[0].reaction_charge_coulombs +
                result.terminals.potential_volts[1] *
                    result.discrete_reactions[1].reaction_charge_coulombs;
            NOVA_REQUIRE_NEAR(terminal_work,
                              2.0 * result.field.stored_energy_joules, 1.0e-8);
        }
    }
}
