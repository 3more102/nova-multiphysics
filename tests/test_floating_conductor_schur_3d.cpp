#include "test_framework.hpp"

#include "nova/physics/floating_conductor_schur_3d.hpp"

#include <algorithm>
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


NOVA_TEST(floating_schur_reference_voltage_shift_and_electrode_permutation) {
    // With only two electrodes and insulated sidewalls, a common potential
    // shift is a gauge transformation: V -> V + constant, but E, D, Q, U
    // and the discrete terminal capacitance are unchanged.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double x, double y, double) {
        return 1.0 + 0.5 * x + 0.125 * y;
    };
    auto groups = plates(problem, box.boundary_faces);
    const auto baseline = solve_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, 5.0}});
    std::reverse(groups.begin(), groups.end());
    const auto shifted = solve_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 4.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, 5.0}});

    NOVA_REQUIRE(baseline.field.linear_solver.converged);
    NOVA_REQUIRE(shifted.field.linear_solver.converged);
    NOVA_REQUIRE(baseline.terminals.electrode_ids ==
                 std::vector<int>({10, 20}));
    NOVA_REQUIRE(shifted.terminals.electrode_ids ==
                 baseline.terminals.electrode_ids);
    NOVA_REQUIRE(baseline.discrete_capacitance.matrix_farads.size() ==
                 shifted.discrete_capacitance.matrix_farads.size());
    for (std::size_t i = 0;
         i < baseline.discrete_capacitance.matrix_farads.size(); ++i) {
        NOVA_REQUIRE_NEAR(baseline.discrete_capacitance.matrix_farads[i],
                          shifted.discrete_capacitance.matrix_farads[i], 1.0e-8);
    }
    NOVA_REQUIRE_NEAR(shifted.terminals.potential_volts[0], 4.0, 1.0e-7);
    NOVA_REQUIRE_NEAR(
        shifted.terminals.potential_volts[1] -
        baseline.terminals.potential_volts[1], 4.0, 2.0e-6);
    NOVA_REQUIRE_NEAR(shifted.field.stored_energy_joules,
                      baseline.field.stored_energy_joules, 2.0e-6);
    NOVA_REQUIRE(baseline.discrete_reactions.size() == 2U);
    NOVA_REQUIRE(shifted.discrete_reactions.size() == 2U);
    for (std::size_t i = 0; i < 2U; ++i) {
        NOVA_REQUIRE_NEAR(
            baseline.discrete_reactions[i].reaction_charge_coulombs,
            shifted.discrete_reactions[i].reaction_charge_coulombs, 2.0e-6);
    }
    NOVA_REQUIRE_NEAR(
        baseline.discrete_reactions[0].reaction_charge_coulombs,
        -5.0, 2.0e-6);
    NOVA_REQUIRE_NEAR(
        shifted.discrete_reactions[1].reaction_charge_coulombs,
        5.0, 2.0e-6);

    NOVA_REQUIRE(baseline.field.potential.size() ==
                 shifted.field.potential.size());
    for (std::size_t i = 0; i < baseline.field.potential.size(); ++i) {
        NOVA_REQUIRE_NEAR(shifted.field.potential[i] -
                          baseline.field.potential[i], 4.0, 2.0e-6);
    }
    NOVA_REQUIRE(baseline.field.electric_field.size() ==
                 shifted.field.electric_field.size());
    NOVA_REQUIRE(baseline.field.displacement.size() ==
                 shifted.field.displacement.size());
    for (std::size_t i = 0; i < baseline.field.electric_field.size(); ++i) {
        const auto& a = baseline.field.electric_field[i];
        const auto& b = shifted.field.electric_field[i];
        const auto& da = baseline.field.displacement[i];
        const auto& db = shifted.field.displacement[i];
        NOVA_REQUIRE_NEAR(a.x, b.x, 2.0e-6);
        NOVA_REQUIRE_NEAR(a.y, b.y, 2.0e-6);
        NOVA_REQUIRE_NEAR(a.z, b.z, 2.0e-6);
        NOVA_REQUIRE_NEAR(da.x, db.x, 2.0e-6);
        NOVA_REQUIRE_NEAR(da.y, db.y, 2.0e-6);
        NOVA_REQUIRE_NEAR(da.z, db.z, 2.0e-6);
    }
}

NOVA_TEST(floating_schur_charge_sign_reverses_fields_not_energy) {
    // Linear source-free electrostatics is odd in applied terminal charge.
    // The quadratic stored field energy must remain unchanged under Q -> -Q.
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double x, double y, double) {
        return 1.0 + 0.5 * x + 0.125 * y;
    };
    const auto groups = plates(problem, box.boundary_faces);
    const auto positive = solve_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, 3.0}});
    const auto negative = solve_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, -3.0}});
    const auto zero = solve_floating_conductor_schur_3d(
        problem, groups, std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, 0.0}});

    NOVA_REQUIRE(positive.field.linear_solver.converged);
    NOVA_REQUIRE(negative.field.linear_solver.converged);
    NOVA_REQUIRE(zero.field.linear_solver.converged);
    NOVA_REQUIRE(positive.field.stored_energy_joules > 0.0);
    NOVA_REQUIRE_NEAR(positive.field.stored_energy_joules,
                      negative.field.stored_energy_joules, 2.0e-6);
    NOVA_REQUIRE_NEAR(zero.field.stored_energy_joules, 0.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(positive.terminals.potential_volts[1],
                      -negative.terminals.potential_volts[1], 2.0e-6);
    NOVA_REQUIRE_NEAR(zero.terminals.potential_volts[1], 0.0, 1.0e-8);
    NOVA_REQUIRE(positive.discrete_reactions.size() == 2U);
    NOVA_REQUIRE(negative.discrete_reactions.size() == 2U);
    for (std::size_t i = 0; i < 2U; ++i) {
        NOVA_REQUIRE_NEAR(
            positive.discrete_reactions[i].reaction_charge_coulombs,
            -negative.discrete_reactions[i].reaction_charge_coulombs, 2.0e-6);
        NOVA_REQUIRE_NEAR(
            zero.discrete_reactions[i].reaction_charge_coulombs,
            0.0, 1.0e-8);
    }
    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(positive.field.potential[i],
                          -negative.field.potential[i], 2.0e-6);
        NOVA_REQUIRE_NEAR(zero.field.potential[i], 0.0, 1.0e-8);
    }
    for (std::size_t i = 0; i < problem.mesh.elements.size(); ++i) {
        const auto& ep = positive.field.electric_field[i];
        const auto& en = negative.field.electric_field[i];
        const auto& dp = positive.field.displacement[i];
        const auto& dn = negative.field.displacement[i];
        NOVA_REQUIRE_NEAR(ep.x, -en.x, 2.0e-6);
        NOVA_REQUIRE_NEAR(ep.y, -en.y, 2.0e-6);
        NOVA_REQUIRE_NEAR(ep.z, -en.z, 2.0e-6);
        NOVA_REQUIRE_NEAR(dp.x, -dn.x, 2.0e-6);
        NOVA_REQUIRE_NEAR(dp.y, -dn.y, 2.0e-6);
        NOVA_REQUIRE_NEAR(dp.z, -dn.z, 2.0e-6);
    }
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
