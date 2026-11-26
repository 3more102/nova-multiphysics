#include "test_framework.hpp"
#include "nova/physics/floating_electrodes_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

using nova::physics::CapacitanceMatrix3D;
using nova::physics::FixedElectrodePotential3D;
using nova::physics::FloatingElectrodeCharge3D;
using nova::physics::solve_floating_electrode_potentials_3d;

template <class Callable>
void require_invalid(const Callable& operation) {
    bool invalid = false;
    try { operation(); }
    catch (const std::invalid_argument&) { invalid = true; }
    NOVA_REQUIRE(invalid);
}

} // namespace

NOVA_TEST(floating_electrodes_parallel_plate_charge_to_voltage) {
    const CapacitanceMatrix3D c{{10, 20}, {2.0, -2.0, -2.0, 2.0}};
    const std::vector<FixedElectrodePotential3D> fixed{{10, 0.0}};
    const std::vector<FloatingElectrodeCharge3D> free{{20, 20.0}};
    const auto result = solve_floating_electrode_potentials_3d(c, fixed, free);
    NOVA_REQUIRE(result.electrode_ids == c.electrode_ids);
    NOVA_REQUIRE(result.potential_volts.size() == 2U);
    NOVA_REQUIRE_NEAR(result.potential_volts[0], 0.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.potential_volts[1], 10.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[0], -20.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[1], 20.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.max_floating_charge_residual_coulombs, 0.0, 1.0e-12);

    const auto shifted = solve_floating_electrode_potentials_3d(
        c, std::vector<FixedElectrodePotential3D>{{20, 4.0}},
        std::vector<FloatingElectrodeCharge3D>{{10, -6.0}});
    NOVA_REQUIRE_NEAR(shifted.potential_volts[0], 1.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(shifted.potential_volts[1], 4.0, 1.0e-12);
}

NOVA_TEST(floating_electrodes_three_terminal_coupling) {
    const CapacitanceMatrix3D c{
        {10, 20, 30},
        {3.0, -1.0, -2.0, -1.0, 2.0, -1.0, -2.0, -1.0, 3.0}
    };
    // Deliberately supply floating terminals in reverse electrode-ID order.
    const auto result = solve_floating_electrode_potentials_3d(
        c, std::vector<FixedElectrodePotential3D>{{10, 3.0}},
        std::vector<FloatingElectrodeCharge3D>{{30, -4.0}, {20, 3.0}});
    NOVA_REQUIRE_NEAR(result.potential_volts[0], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.potential_volts[1], 4.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.potential_volts[2], 2.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[0], 1.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[1], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[2], -4.0, 1.0e-12);
}

NOVA_TEST(floating_electrodes_external_ground_and_all_fixed) {
    // A positive-definite terminal matrix models an independent grounded
    // reference; unlike an isolated Maxwell matrix, all terminals may float.
    const CapacitanceMatrix3D c{{1, 2}, {3.0, -1.0, -1.0, 2.0}};
    const auto all_free = solve_floating_electrode_potentials_3d(
        c, std::vector<FixedElectrodePotential3D>{},
        std::vector<FloatingElectrodeCharge3D>{{1, 3.0}, {2, 4.0}});
    NOVA_REQUIRE_NEAR(all_free.potential_volts[0], 2.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_free.potential_volts[1], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_free.estimated_charge_coulombs[0], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_free.estimated_charge_coulombs[1], 4.0, 1.0e-12);

    const auto all_fixed = solve_floating_electrode_potentials_3d(
        c, std::vector<FixedElectrodePotential3D>{{2, 3.0}, {1, 2.0}},
        std::vector<FloatingElectrodeCharge3D>{});
    NOVA_REQUIRE_NEAR(all_fixed.estimated_charge_coulombs[0], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_fixed.estimated_charge_coulombs[1], 4.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_fixed.max_floating_charge_residual_coulombs, 0.0, 1.0e-12);
}

NOVA_TEST(floating_electrodes_rejects_invalid_inputs) {
    const CapacitanceMatrix3D c{{10, 20}, {2.0, -2.0, -2.0, 2.0}};
    const std::vector<FixedElectrodePotential3D> fixed{{10, 0.0}};
    const std::vector<FloatingElectrodeCharge3D> free{{20, 20.0}};
    const auto run = [&](const CapacitanceMatrix3D& matrix,
                         const std::vector<FixedElectrodePotential3D>& k,
                         const std::vector<FloatingElectrodeCharge3D>& f) {
        (void)solve_floating_electrode_potentials_3d(matrix, k, f);
    };

    require_invalid([&] { run({{}, {}}, fixed, free); });
    require_invalid([&] { run({{10, 20}, {2.0}}, fixed, free); });
    require_invalid([&] { run({{20, 10}, c.matrix_farads}, fixed, free); });
    require_invalid([&] { run({{10, 10}, c.matrix_farads}, fixed, free); });
    require_invalid([&] { run({{-1, 20}, c.matrix_farads}, fixed, free); });
    require_invalid([&] {
        run({{10, 20}, {2.0, -2.0, -2.0, std::numeric_limits<double>::infinity()}},
            fixed, free);
    });
    require_invalid([&] {
        run(c, {{10, std::numeric_limits<double>::quiet_NaN()}}, free);
    });
    require_invalid([&] {
        run(c, fixed, {{20, std::numeric_limits<double>::infinity()}});
    });
    require_invalid([&] { run(c, {{999, 0.0}}, free); });
    require_invalid([&] { run(c, {{10, 0.0}, {10, 0.0}}, free); });
    require_invalid([&] { run(c, fixed, {{20, 1.0}, {20, 1.0}}); });
    require_invalid([&] { run(c, fixed, {{10, 0.0}}); });
    require_invalid([&] { run(c, fixed, {}); });

    // Isolated two-electrode Maxwell matrix has a constant-voltage nullspace.
    require_invalid([&] {
        run(c, {}, {{10, -20.0}, {20, 20.0}});
    });
    // A badly conditioned floating block is not silently inverted.
    require_invalid([&] {
        run({{10, 20}, {1.0, 0.0, 0.0, 1.0e-14}}, {},
            {{10, 1.0}, {20, 1.0}});
    });
}


NOVA_TEST(affine_floating_electrodes_nonzero_offset_and_source_free_equivalence) {
    const CapacitanceMatrix3D c{{10, 20}, {2.0, -2.0, -2.0, 2.0}};
    const std::vector<FixedElectrodePotential3D> fixed{{10, 0.0}};
    const std::vector<FloatingElectrodeCharge3D> loaded{{20, 23.0}};
    const auto result = nova::physics::solve_affine_floating_electrode_potentials_3d(
        c, std::vector<double>{-1.0, 3.0}, fixed, loaded);
    // Q_high=2*V_high + 3 C, so Q_high=23 C implies V_high=10 V.
    NOVA_REQUIRE(result.electrode_ids == c.electrode_ids);
    NOVA_REQUIRE_NEAR(result.potential_volts[0], 0.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.potential_volts[1], 10.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[0], -21.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.estimated_charge_coulombs[1], 23.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(result.max_floating_charge_residual_coulombs, 0.0, 1.0e-12);

    const auto source_free = solve_floating_electrode_potentials_3d(
        c, fixed, std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    const auto zero_offset = nova::physics::solve_affine_floating_electrode_potentials_3d(
        c, std::vector<double>{0.0, 0.0}, fixed,
        std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    NOVA_REQUIRE(zero_offset.potential_volts == source_free.potential_volts);
    NOVA_REQUIRE(zero_offset.estimated_charge_coulombs ==
                 source_free.estimated_charge_coulombs);
}

NOVA_TEST(affine_floating_electrodes_coupled_three_terminal_and_all_fixed) {
    const CapacitanceMatrix3D c{
        {10, 20, 30},
        {3.0, -1.0, -2.0, -1.0, 2.0, -1.0, -2.0, -1.0, 3.0}
    };
    // V=(3,4,2) gives C*V=(1,3,-4), Q_zero=(2,-1,3),
    // and Q_total=(3,2,-1). Floating requests are in reverse ID order.
    const auto coupled = nova::physics::solve_affine_floating_electrode_potentials_3d(
        c, std::vector<double>{2.0, -1.0, 3.0},
        std::vector<FixedElectrodePotential3D>{{10, 3.0}},
        std::vector<FloatingElectrodeCharge3D>{{30, -1.0}, {20, 2.0}});
    NOVA_REQUIRE_NEAR(coupled.potential_volts[0], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(coupled.potential_volts[1], 4.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(coupled.potential_volts[2], 2.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(coupled.estimated_charge_coulombs[0], 3.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(coupled.estimated_charge_coulombs[1], 2.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(coupled.estimated_charge_coulombs[2], -1.0, 1.0e-12);

    const auto all_fixed = nova::physics::solve_affine_floating_electrode_potentials_3d(
        CapacitanceMatrix3D{{1, 2}, {3.0, -1.0, -1.0, 2.0}},
        std::vector<double>{-1.0, 3.0},
        std::vector<FixedElectrodePotential3D>{{2, 3.0}, {1, 2.0}},
        std::vector<FloatingElectrodeCharge3D>{});
    NOVA_REQUIRE_NEAR(all_fixed.estimated_charge_coulombs[0], 2.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_fixed.estimated_charge_coulombs[1], 7.0, 1.0e-12);
    NOVA_REQUIRE_NEAR(all_fixed.max_floating_charge_residual_coulombs, 0.0, 1.0e-12);
}

NOVA_TEST(affine_floating_electrodes_rejects_invalid_offsets_and_constraints) {
    const CapacitanceMatrix3D c{{10, 20}, {2.0, -2.0, -2.0, 2.0}};
    const std::vector<FixedElectrodePotential3D> fixed{{10, 0.0}};
    const std::vector<FloatingElectrodeCharge3D> floating{{20, 20.0}};
    const auto solve = [&](const CapacitanceMatrix3D& matrix,
                           const std::vector<double>& offsets,
                           const std::vector<FixedElectrodePotential3D>& f,
                           const std::vector<FloatingElectrodeCharge3D>& q) {
        (void)nova::physics::solve_affine_floating_electrode_potentials_3d(
            matrix, offsets, f, q);
    };
    require_invalid([&] { solve(c, {0.0}, fixed, floating); });
    require_invalid([&] { solve(c, {0.0, 0.0, 0.0}, fixed, floating); });
    require_invalid([&] {
        solve(c, {std::numeric_limits<double>::quiet_NaN(), 0.0}, fixed, floating);
    });
    require_invalid([&] {
        solve(c, {0.0, std::numeric_limits<double>::infinity()}, fixed, floating);
    });
    require_invalid([&] { solve(c, {0.0, 0.0}, fixed, {{99, 1.0}}); });
    require_invalid([&] {
        solve(c, {0.0, 0.0}, fixed, {{20, 1.0}, {20, 1.0}});
    });
    require_invalid([&] {
        solve(c, {0.0, 0.0}, {}, {{10, -20.0}, {20, 20.0}});
    });
}

NOVA_TEST(floating_electrodes_reconstructs_3d_field_from_charge) {
    using nova::physics::Electrostatics3DProblem;
    using nova::physics::ElectrodeBoundaryGroup3D;
    using nova::physics::reconstruct_floating_electrostatic_field_3d;
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    std::vector<ElectrodeBoundaryGroup3D> groups{{10, {}}, {20, {}}};
    for (const auto& face : box.boundary_faces) {
        const double x = problem.mesh.nodes[face[0]].x;
        if (problem.mesh.nodes[face[1]].x == x &&
            problem.mesh.nodes[face[2]].x == x) {
            if (x == 0.0) groups[0].faces.push_back(face);
            if (x == 1.0) {
                // Ignore caller face winding when evaluating the final Q.
                groups[1].faces.push_back({face[2], face[1], face[0]});
            }
        }
    }
    NOVA_REQUIRE(!groups[0].faces.empty() && !groups[1].faces.empty());
    const auto result = reconstruct_floating_electrostatic_field_3d(
        problem, groups,
        std::vector<FixedElectrodePotential3D>{{10, 0.0}},
        std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    NOVA_REQUIRE(result.field.linear_solver.converged);
    NOVA_REQUIRE(result.capacitance.electrode_ids.size() == 2U);
    NOVA_REQUIRE_NEAR(result.terminals.potential_volts[1], 10.0, 1.0e-8);
    NOVA_REQUIRE(result.electrode_charges.size() == 2U);
    NOVA_REQUIRE_NEAR(
        result.electrode_charges[0].conductor_charge_coulombs, -20.0, 1.0e-8);
    NOVA_REQUIRE_NEAR(
        result.electrode_charges[1].conductor_charge_coulombs, 20.0, 1.0e-8);
    NOVA_REQUIRE(result.max_floating_charge_mismatch_coulombs < 1.0e-8);
    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(
            result.field.potential[i], 10.0 * problem.mesh.nodes[i].x, 1.0e-8);
    }
    for (const auto& field : result.field.electric_field) {
        NOVA_REQUIRE_NEAR(field.x, -10.0, 1.0e-8);
        NOVA_REQUIRE_NEAR(field.y, 0.0, 1.0e-8);
        NOVA_REQUIRE_NEAR(field.z, 0.0, 1.0e-8);
    }
    NOVA_REQUIRE_NEAR(result.field.stored_energy_joules, 100.0, 1.0e-7);
}

NOVA_TEST(floating_electrodes_reconstructs_field_with_external_ground) {
    using nova::physics::Electrostatics3DProblem;
    using nova::physics::ElectrodeBoundaryGroup3D;
    using nova::physics::reconstruct_floating_electrostatic_field_3d;
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    ElectrodeBoundaryGroup3D floating{20, {}};
    for (const auto& face : box.boundary_faces) {
        const double x = problem.mesh.nodes[face[0]].x;
        if (problem.mesh.nodes[face[1]].x != x ||
            problem.mesh.nodes[face[2]].x != x) continue;
        if (x == 0.0) {
            // Capacitance extraction treats every external reference as 0 V,
            // regardless of the originally supplied Dirichlet callback.
            problem.potential_faces.push_back({
                face, [](double, double, double) { return 7.0; }});
        }
        if (x == 1.0) floating.faces.push_back(face);
    }
    const auto result = reconstruct_floating_electrostatic_field_3d(
        problem, std::vector<ElectrodeBoundaryGroup3D>{floating},
        std::vector<FixedElectrodePotential3D>{},
        std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    NOVA_REQUIRE(result.field.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.terminals.potential_volts[0], 10.0, 1.0e-8);
    NOVA_REQUIRE(result.max_floating_charge_mismatch_coulombs < 1.0e-8);
    for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
        NOVA_REQUIRE_NEAR(
            result.field.potential[i], 10.0 * problem.mesh.nodes[i].x, 1.0e-8);
    }
    NOVA_REQUIRE_NEAR(result.field.stored_energy_joules, 100.0, 1.0e-7);
}

NOVA_TEST(floating_electrodes_reconstruction_rejects_impressed_loads) {
    using nova::physics::Electrostatics3DProblem;
    using nova::physics::ElectrodeBoundaryGroup3D;
    using nova::physics::reconstruct_floating_electrostatic_field_3d;
    const auto box = nova::fem::make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1U, 1U, 1U);
    Electrostatics3DProblem problem;
    problem.mesh = box.mesh;
    problem.permittivity = [](double, double, double) { return 2.0; };
    std::vector<ElectrodeBoundaryGroup3D> groups{{10, {}}, {20, {}}};
    for (const auto& face : box.boundary_faces) {
        const double x = problem.mesh.nodes[face[0]].x;
        if (problem.mesh.nodes[face[1]].x != x ||
            problem.mesh.nodes[face[2]].x != x) continue;
        if (x == 0.0) groups[0].faces.push_back(face);
        if (x == 1.0) groups[1].faces.push_back(face);
    }
    const auto run = [&] {
        (void)reconstruct_floating_electrostatic_field_3d(
            problem, groups,
            std::vector<FixedElectrodePotential3D>{{10, 0.0}},
            std::vector<FloatingElectrodeCharge3D>{{20, 20.0}});
    };
    problem.free_charge_density = [](double, double, double) { return 1.0; };
    require_invalid(run);
    problem.free_charge_density = {};
    problem.cell_free_charge_density.assign(problem.mesh.elements.size(), 1.0);
    require_invalid(run);
    problem.cell_free_charge_density.clear();
    problem.displacement_flux_faces.push_back({
        box.boundary_faces.front(),
        [](double, double, double) { return 1.0; }});
    require_invalid(run);
}
