#include "test_framework.hpp"

#include "nova/physics/floating_conductor_schur_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

// Analytical manufactured capacitor with epsilon(x) = 1+x F/m,
// unit cross-sectional area, x in [0,1], left plate V=0 and right plate
// charge Q=+1 C. The unique continuum solution is
//   V(x) = log(1+x) V,  V_right = log(2) V,
//   E_x = -1/(1+x) V/m,  U = log(2)/2 J,
//   C_equivalent = 1/log(2) F.
// The side walls are naturally insulated; no volume charge is present.
NOVA_TEST(floating_schur_dielectric_capacitor_mesh_convergence) {
    const double exact_voltage = std::log(2.0);
    const double exact_capacitance = 1.0 / exact_voltage;
    const double exact_energy = 0.5 * exact_voltage;

    const std::array<std::size_t, 3> resolutions{1U, 2U, 4U};
    std::array<double, 3> voltage_errors{};
    std::array<double, 3> capacitance_errors{};
    std::array<double, 3> field_node_errors{};

    for (std::size_t level = 0; level < resolutions.size(); ++level) {
        const auto n = resolutions[level];
        const auto box = nova::fem::make_box_mesh_3d(
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, n, n);
        nova::physics::Electrostatics3DProblem problem;
        problem.mesh = box.mesh;
        problem.permittivity = [](double x, double, double) {
            return 1.0 + x;
        };
        std::vector<nova::physics::ElectrodeBoundaryGroup3D> groups{
            {10, {}}, {20, {}}
        };
        for (const auto& face : box.boundary_faces) {
            const double x0 = problem.mesh.nodes[face[0]].x;
            const double x1 = problem.mesh.nodes[face[1]].x;
            const double x2 = problem.mesh.nodes[face[2]].x;
            if (x0 == 0.0 && x1 == 0.0 && x2 == 0.0) {
                groups[0].faces.push_back(face);
            } else if (x0 == 1.0 && x1 == 1.0 && x2 == 1.0) {
                groups[1].faces.push_back(face);
            }
        }
        NOVA_REQUIRE(!groups[0].faces.empty() && !groups[1].faces.empty());

        const auto result = nova::physics::solve_floating_conductor_schur_3d(
            problem, groups,
            std::vector<nova::physics::FixedElectrodePotential3D>{{10, 0.0}},
            std::vector<nova::physics::FloatingElectrodeCharge3D>{{20, 1.0}});
        NOVA_REQUIRE(result.field.linear_solver.converged);
        NOVA_REQUIRE(result.discrete_capacitance.matrix_farads.size() == 4U);
        NOVA_REQUIRE(result.discrete_reactions.size() == 2U);
        NOVA_REQUIRE_NEAR(
            result.discrete_reactions[1].reaction_charge_coulombs, 1.0, 1.0e-7);
        NOVA_REQUIRE_NEAR(
            result.discrete_reactions[0].reaction_charge_coulombs, -1.0, 1.0e-7);
        NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-7);
        NOVA_REQUIRE_NEAR(
            result.terminals.potential_volts[1] *
                result.discrete_reactions[1].reaction_charge_coulombs,
            2.0 * result.field.stored_energy_joules, 1.0e-6);

        const double numerical_voltage = result.terminals.potential_volts[1];
        voltage_errors[level] = std::abs(numerical_voltage - exact_voltage);
        capacitance_errors[level] = std::abs(
            result.discrete_capacitance.matrix_farads[3] - exact_capacitance);
        double maximum_nodal_error = 0.0;
        for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
            const double analytic = std::log1p(problem.mesh.nodes[i].x);
            maximum_nodal_error = std::max(
                maximum_nodal_error,
                std::abs(result.field.potential[i] - analytic));
        }
        field_node_errors[level] = maximum_nodal_error;
        if (level + 1U == resolutions.size()) {
            NOVA_REQUIRE_NEAR(result.field.stored_energy_joules,
                              exact_energy, 0.01);
        }
    }
    // Continuum convergence is a distinct check from satisfying the FEM
    // charge constraints. These are errors against an analytic solution.
    NOVA_REQUIRE(voltage_errors[2] < voltage_errors[1]);
    NOVA_REQUIRE(voltage_errors[1] < voltage_errors[0]);
    NOVA_REQUIRE(capacitance_errors[2] < capacitance_errors[1]);
    NOVA_REQUIRE(capacitance_errors[1] < capacitance_errors[0]);
    NOVA_REQUIRE(field_node_errors[2] < field_node_errors[0]);
    NOVA_REQUIRE(voltage_errors[2] < 0.01);
}


/*
 * Independently soluble SERIES dielectric capacitor with an interface
 * conforming to the tetrahedral mesh: epsilon = 2 F/m for x < 1/2,
 * epsilon = 6 F/m for x > 1/2, cross-section A = 1 m^2.
 *
 * Ground x=0 and prescribe Q=+1 C at the floating x=1 electrode.
 * Gauss's law gives constant D_x=-1 C/m^2 (and insulated sidewalls).
 * V(x) = x/2 for x<=1/2; V(x)=1/4+(x-1/2)/6 otherwise.
 * V_high=1/3 V, C=3 F, U=1/6 J. The interface jump in E is
 * physical, while D normal is continuous. The exact potential is
 * representable by conforming P1 tetrahedra at both mesh resolutions.
 */
NOVA_TEST(floating_schur_two_layer_dielectric_interface_exact_solution) {
    for (const std::size_t n : {2U, 4U}) {
        const auto box = nova::fem::make_box_mesh_3d(
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, n, n);
        nova::physics::Electrostatics3DProblem problem;
        problem.mesh = box.mesh;
        problem.permittivity = [](double x, double, double) {
            return x < 0.5 ? 2.0 : 6.0;
        };

        std::vector<nova::physics::ElectrodeBoundaryGroup3D> electrodes{
            {10, {}}, {20, {}}
        };
        for (const auto& face : box.boundary_faces) {
            const auto x0 = problem.mesh.nodes[face[0]].x;
            const auto x1 = problem.mesh.nodes[face[1]].x;
            const auto x2 = problem.mesh.nodes[face[2]].x;
            if (x0 == 0.0 && x1 == 0.0 && x2 == 0.0) {
                electrodes[0].faces.push_back(face);
            } else if (x0 == 1.0 && x1 == 1.0 && x2 == 1.0) {
                // Check that electrode charge is independent of face winding.
                electrodes[1].faces.push_back({face[2], face[1], face[0]});
            }
        }
        NOVA_REQUIRE(!electrodes[0].faces.empty());
        NOVA_REQUIRE(!electrodes[1].faces.empty());

        const auto result = nova::physics::solve_floating_conductor_schur_3d(
            problem, electrodes,
            std::vector<nova::physics::FixedElectrodePotential3D>{{10, 0.0}},
            std::vector<nova::physics::FloatingElectrodeCharge3D>{{20, 1.0}});
        NOVA_REQUIRE(result.field.linear_solver.converged);
        NOVA_REQUIRE(result.discrete_capacitance.electrode_ids ==
                     std::vector<int>({10, 20}));
        const auto& capacitance = result.discrete_capacitance.matrix_farads;
        NOVA_REQUIRE(capacitance.size() == 4U);
        NOVA_REQUIRE_NEAR(capacitance[0], 3.0, 1.0e-7);
        NOVA_REQUIRE_NEAR(capacitance[1], -3.0, 1.0e-7);
        NOVA_REQUIRE_NEAR(capacitance[2], -3.0, 1.0e-7);
        NOVA_REQUIRE_NEAR(capacitance[3], 3.0, 1.0e-7);
        NOVA_REQUIRE(result.discrete_capacitance.max_reciprocity_defect_farads <
                     1.0e-7);

        NOVA_REQUIRE_NEAR(result.terminals.potential_volts[0], 0.0, 1.0e-8);
        NOVA_REQUIRE_NEAR(result.terminals.potential_volts[1], 1.0 / 3.0, 1.0e-7);
        NOVA_REQUIRE(result.discrete_reactions.size() == 2U);
        NOVA_REQUIRE_NEAR(result.discrete_reactions[0].reaction_charge_coulombs,
                          -1.0, 1.0e-7);
        NOVA_REQUIRE_NEAR(result.discrete_reactions[1].reaction_charge_coulombs,
                          1.0, 1.0e-7);
        NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-7);
        NOVA_REQUIRE(result.max_flux_reaction_disagreement_coulombs < 1.0e-7);
        NOVA_REQUIRE_NEAR(result.field.stored_energy_joules, 1.0 / 6.0, 1.0e-7);

        for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
            const double x = problem.mesh.nodes[i].x;
            const double exact = x <= 0.5
                ? x / 2.0 : 0.25 + (x - 0.5) / 6.0;
            NOVA_REQUIRE_NEAR(result.field.potential[i], exact, 1.0e-7);
        }
        NOVA_REQUIRE(result.field.displacement.size() ==
                     problem.mesh.elements.size());
        NOVA_REQUIRE(result.field.electric_field.size() ==
                     problem.mesh.elements.size());
        for (std::size_t i = 0; i < problem.mesh.elements.size(); ++i) {
            const auto& cell = problem.mesh.elements[i];
            const double cx =
                (problem.mesh.nodes[cell[0]].x + problem.mesh.nodes[cell[1]].x +
                 problem.mesh.nodes[cell[2]].x + problem.mesh.nodes[cell[3]].x) /
                4.0;
            const double epsilon = cx < 0.5 ? 2.0 : 6.0;
            NOVA_REQUIRE_NEAR(result.field.electric_field[i].x,
                              -1.0 / epsilon, 1.0e-7);
            NOVA_REQUIRE_NEAR(result.field.displacement[i].x, -1.0, 1.0e-7);
            NOVA_REQUIRE_NEAR(result.field.displacement[i].y, 0.0, 1.0e-7);
            NOVA_REQUIRE_NEAR(result.field.displacement[i].z, 0.0, 1.0e-7);
        }
    }
}

} // namespace
