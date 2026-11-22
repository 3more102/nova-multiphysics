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


// Exact transmission problem with a discontinuous dielectric coefficient.
// epsilon = 2 F/m on 0 <= x < 1/2 and epsilon = 4 F/m on 1/2 < x <= 1.
// The mesh contains the interface plane, so piecewise affine V belongs to P1.
// With unit cross-section, Q_right=+1 C, and V_left=0, the displacement
// D_x=-1 C/m^2 is continuous while E_x jumps from -1/2 to -1/4 V/m.
//   V(x) = x/2                             (x <= 1/2)
//          1/4 + (x-1/2)/4                 (x >= 1/2)
//   V_right = 3/8 V, C_equiv = 8/3 F, U = 3/16 J.
// Run with independently populated cell materials and the spatial callback
// to catch material-mapping errors at internal interfaces.
NOVA_TEST(floating_schur_two_layer_dielectric_exact_interface_patch) {
    constexpr double exact_voltage = 3.0 / 8.0;
    constexpr double exact_capacitance = 8.0 / 3.0;
    constexpr double exact_energy = 3.0 / 16.0;
    const std::array<std::size_t, 2> resolutions{2U, 4U};

    for (const auto n : resolutions) {
        const auto box = nova::fem::make_box_mesh_3d(
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, n, n);
        std::vector<nova::physics::ElectrodeBoundaryGroup3D> groups{
            {10, {}}, {20, {}}
        };
        for (const auto& face : box.boundary_faces) {
            bool left = true;
            bool right = true;
            for (const auto node : face) {
                left = left && box.mesh.nodes[node].x == 0.0;
                right = right && box.mesh.nodes[node].x == 1.0;
            }
            if (left) {
                groups[0].faces.push_back(face);
            } else if (right) {
                groups[1].faces.push_back(face);
            }
        }
        NOVA_REQUIRE(!groups[0].faces.empty() && !groups[1].faces.empty());

        nova::physics::Electrostatics3DProblem problem;
        problem.mesh = box.mesh;
        problem.cell_permittivity.reserve(problem.mesh.elements.size());
        for (const auto& tetra : problem.mesh.elements) {
            double centroid_x = 0.0;
            for (const auto node : tetra) {
                centroid_x += problem.mesh.nodes[node].x / 4.0;
            }
            NOVA_REQUIRE(std::abs(centroid_x - 0.5) > 1.0e-12);
            problem.cell_permittivity.push_back(centroid_x < 0.5 ? 2.0 : 4.0);
        }

        for (int material_mode = 0; material_mode < 2; ++material_mode) {
            if (material_mode == 1) {
                problem.cell_permittivity.clear();
                problem.permittivity = [](double x, double, double) {
                    return x < 0.5 ? 2.0 : 4.0;
                };
            }

            const auto result = nova::physics::solve_floating_conductor_schur_3d(
                problem, groups,
                std::vector<nova::physics::FixedElectrodePotential3D>{{10, 0.0}},
                std::vector<nova::physics::FloatingElectrodeCharge3D>{{20, 1.0}});
            NOVA_REQUIRE(result.field.linear_solver.converged);
            NOVA_REQUIRE(result.discrete_capacitance.matrix_farads.size() == 4U);
            NOVA_REQUIRE_NEAR(
                result.discrete_capacitance.matrix_farads[3],
                exact_capacitance, 1.0e-6);
            NOVA_REQUIRE_NEAR(
                result.terminals.potential_volts[1], exact_voltage, 1.0e-7);
            NOVA_REQUIRE_NEAR(
                result.field.stored_energy_joules, exact_energy, 1.0e-7);
            NOVA_REQUIRE_NEAR(
                result.discrete_reactions[0].reaction_charge_coulombs,
                -1.0, 1.0e-7);
            NOVA_REQUIRE_NEAR(
                result.discrete_reactions[1].reaction_charge_coulombs,
                1.0, 1.0e-7);
            NOVA_REQUIRE(result.max_floating_reaction_mismatch_coulombs < 1.0e-7);
            NOVA_REQUIRE(result.max_flux_reaction_disagreement_coulombs < 1.0e-7);

            for (std::size_t i = 0; i < problem.mesh.nodes.size(); ++i) {
                const double x = problem.mesh.nodes[i].x;
                const double expected = x <= 0.5
                    ? 0.5 * x
                    : 0.25 + 0.25 * (x - 0.5);
                NOVA_REQUIRE_NEAR(result.field.potential[i], expected, 1.0e-7);
            }
            NOVA_REQUIRE(
                result.field.electric_field.size() ==
                problem.mesh.elements.size());
            NOVA_REQUIRE(
                result.field.displacement.size() ==
                problem.mesh.elements.size());
            for (std::size_t i = 0; i < problem.mesh.elements.size(); ++i) {
                const auto& tetra = problem.mesh.elements[i];
                double centroid_x = 0.0;
                for (const auto node : tetra) {
                    centroid_x += problem.mesh.nodes[node].x / 4.0;
                }
                const double expected_ex = centroid_x < 0.5 ? -0.5 : -0.25;
                const auto& electric = result.field.electric_field[i];
                const auto& displacement = result.field.displacement[i];
                NOVA_REQUIRE_NEAR(electric.x, expected_ex, 1.0e-7);
                NOVA_REQUIRE_NEAR(electric.y, 0.0, 1.0e-7);
                NOVA_REQUIRE_NEAR(electric.z, 0.0, 1.0e-7);
                NOVA_REQUIRE_NEAR(displacement.x, -1.0, 1.0e-7);
                NOVA_REQUIRE_NEAR(displacement.y, 0.0, 1.0e-7);
                NOVA_REQUIRE_NEAR(displacement.z, 0.0, 1.0e-7);
            }
            NOVA_REQUIRE_NEAR(
                result.terminals.potential_volts[1] *
                    result.discrete_reactions[1].reaction_charge_coulombs,
                2.0 * result.field.stored_energy_joules, 1.0e-7);
        }
    }
}

} // namespace
