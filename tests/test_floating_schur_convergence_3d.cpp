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

} // namespace
