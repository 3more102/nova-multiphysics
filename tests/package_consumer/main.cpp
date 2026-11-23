#include "nova/physics/electrothermal_3d.hpp"
#include "nova/physics/nonlinear_electrothermal_3d.hpp"
#include "nova/physics/transient_electrothermal_3d.hpp"
#include "nova/physics/adaptive_transient_electrothermal_3d.hpp"
#include "nova/fem/adaptive_diffusion_3d.hpp"
#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/poisson_1d.hpp"
#include "nova/fem/residual_error_estimator_3d.hpp"
#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/conjugate_gradient.hpp"
#include "nova/linalg/incomplete_cholesky.hpp"
#include "nova/physics/electrostatics_3d.hpp"
#include "nova/physics/electrode_charge_3d.hpp"
#include "nova/physics/floating_electrodes_3d.hpp"
#include "nova/physics/floating_conductor_schur_3d.hpp"
#include "nova/io/gmsh_tetrahedral.hpp"
#include "nova/io/thermal_checkpoint_3d.hpp"

#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

int main() {
    // Exercise 1D FEM and CSR algebra across the installed library boundary.
    nova::fem::Poisson1DProblem p;
    p.elements = 16;
    p.source = [](double) { return 2.0; };
    const auto one_d = nova::fem::solve_poisson_1d(p);
    if (one_d.values.size() != 17 ||
        std::abs(one_d.values[8] - 0.25) > 1.0e-7) return 1;

    const auto a = nova::linalg::CsrMatrix::from_triplets(2, 2,
        {{0, 0, 2.0}, {1, 1, 3.0}});
    const std::vector<double> x{4.0, 5.0};
    const auto y = a.multiply(x);
    if (y.size() != 2 || std::abs(y[0] - 8.0) > 1.0e-12 ||
        std::abs(y[1] - 15.0) > 1.0e-12) return 2;

    // Exercise the installed IC(0) factor and CG selection API.
    const auto spd = nova::linalg::CsrMatrix::from_triplets(2, 2,
        {{0, 0, 4.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 3.0}});
    const std::array<double, 2> rhs{1.0, 2.0};
    nova::linalg::SolverOptions ic_options;
    ic_options.preconditioner = nova::linalg::PreconditionerType::IncompleteCholesky0;
    ic_options.relative_tolerance = 1.0e-12;
    ic_options.absolute_tolerance = 1.0e-14;
    const auto ic_solution =
        nova::linalg::solve_conjugate_gradient(spd, rhs, ic_options);
    if (!ic_solution.converged || ic_solution.iterations != 1 ||
        std::abs(ic_solution.solution[0] - 1.0 / 11.0) > 1.0e-12 ||
        std::abs(ic_solution.solution[1] - 7.0 / 11.0) > 1.0e-12) return 8;

    nova::linalg::CachedConjugateGradient cached_cg;
    const std::array<double, 2> second_exact{-1.0, 2.0};
    const auto cached_first = cached_cg.solve(spd, rhs, ic_options);
    const auto cached_second =
        cached_cg.solve(spd, spd.multiply(second_exact), ic_options);
    if (!cached_first.converged || !cached_second.converged ||
        cached_cg.factor_builds() != 1 ||
        std::abs(cached_second.solution[0] + 1.0) > 1.0e-12 ||
        std::abs(cached_second.solution[1] - 2.0) > 1.0e-12) return 9;

    // Independent installed 3D API consumer: -Delta(u)=0, u=x, with
    // fixed temperatures at x=0 and x=1 and insulated remaining faces.
    using namespace nova::fem;
    const auto box = make_box_mesh_3d(
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2U, 2U, 2U);
    ScalarDiffusion3DProblem problem;
    problem.mesh = box.mesh;
    problem.diffusivity = [](double, double, double) { return 1.0; };
    problem.source = [](double, double, double) { return 0.0; };
    for (const auto& face : box.boundary_faces) {
        const bool left = problem.mesh.nodes[face[0]].x == 0.0 &&
            problem.mesh.nodes[face[1]].x == 0.0 &&
            problem.mesh.nodes[face[2]].x == 0.0;
        const bool right = problem.mesh.nodes[face[0]].x == 1.0 &&
            problem.mesh.nodes[face[1]].x == 1.0 &&
            problem.mesh.nodes[face[2]].x == 1.0;
        if (left || right) {
            problem.dirichlet_faces.push_back({face,
                [](double px, double, double) { return px; }});
        }
    }
    const auto three_d = solve_scalar_diffusion_3d(problem);
    if (!three_d.linear_solver.converged ||
        three_d.values.size() != problem.mesh.nodes.size()) return 3;
    for (std::size_t i = 0; i < three_d.values.size(); ++i) {
        if (std::abs(three_d.values[i] - problem.mesh.nodes[i].x) > 1.0e-9)
            return 4;
    }
    // Exercise the installed stationary one-way electrothermal public API.
    nova::physics::Electrothermal3DProblem electrothermal;
    electrothermal.electrical = problem;
    electrothermal.thermal = problem;
    electrothermal.thermal.diffusivity =
        [](double, double, double) { return 4.0; };
    electrothermal.thermal.source = {};
    for (auto& fixed : electrothermal.thermal.dirichlet_faces) {
        fixed.value = [](double, double, double) { return 300.0; };
    }
    const auto coupled = nova::physics::solve_one_way_electrothermal_3d(electrothermal);
    if (!coupled.electrical_solution.linear_solver.converged ||
        !coupled.thermal_solution.linear_solver.converged ||
        std::abs(coupled.total_joule_power_watts - 1.0) > 1.0e-8 ||
        std::abs(coupled.thermal_balance.generated_heat - 1.0) > 1.0e-8 ||
        std::abs(coupled.thermal_balance.signed_balance_defect) > 1.0e-7)
        return 23;

    // Exercise the installed nonlinear electrothermal symbol and temperature
    // feedback independently of the public one-way solve.
    nova::physics::NonlinearElectrothermal3DProblem feedback;
    feedback.base = electrothermal;
    feedback.base.electrical.diffusivity = {};
    feedback.base.electrical.cell_diffusivity.assign(
        box.mesh.elements.size(), 1.0);
    feedback.conductivity_scale = [](double) { return 1.0; };
    const auto nonlinear = nova::physics::solve_nonlinear_electrothermal_3d(feedback);
    if (!nonlinear.converged ||
        nonlinear.iterations < 2 ||
        std::abs(nonlinear.fields.total_joule_power_watts - 1.0) > 1.0e-8 ||
        std::abs(nonlinear.fields.thermal_balance.signed_balance_defect) > 1.0e-7)
        return 24;

    // Optional temperature-dependent k(T) must remain usable from the
    // installed public C++20 SDK, not only from in-tree tests.
    feedback.thermal_conductivity_scale = [](double t) {
        return 1.0 + 0.01 * (t - 300.0);
    };
    // k(T) intentionally requires complete cellwise reference conductivity:
    // the prior SDK exercise used a spatial callback for constant k.
    feedback.base.thermal.diffusivity = {};
    feedback.base.thermal.cell_diffusivity.assign(
        box.mesh.elements.size(), 4.0);
    const auto dual_material =
        nova::physics::solve_nonlinear_electrothermal_3d(feedback);
    if (!dual_material.converged ||
        dual_material.fields.heated_thermal_problem.cell_diffusivity.size() !=
            box.mesh.elements.size() ||
        std::abs(dual_material.fields.total_joule_power_watts - 1.0) > 1.0e-8 ||
        std::abs(dual_material.fields.thermal_balance.signed_balance_defect) > 1.0e-7)
        return 25;

    // Exercise optional adaptive relaxation and public per-iteration history
    // through the installed SDK, including diagnostic field consistency.
    nova::physics::NonlinearElectrothermal3DOptions adaptive_options;
    adaptive_options.adaptive_relaxation = true;
    adaptive_options.minimum_relaxation = 0.1;
    const auto diagnosed = nova::physics::solve_nonlinear_electrothermal_3d(
        feedback, adaptive_options);
    if (!diagnosed.converged ||
        diagnosed.history.size() != diagnosed.iterations ||
        diagnosed.history.size() < 2 ||
        diagnosed.history.back().relaxation_applied != 0.0 ||
        std::abs(diagnosed.history.back().max_normalized_residual -
                 diagnosed.max_normalized_residual) > 1.0e-15)
        return 26;

    // Installed C++20 SDK: transient electrothermal joule source + BE heat.
    nova::physics::TransientElectrothermal3DProblem transient_coupling;
    transient_coupling.electrical = electrothermal.electrical;
    transient_coupling.thermal.diffusion = electrothermal.thermal;
    transient_coupling.thermal.cell_volumetric_capacity.assign(
        box.mesh.elements.size(), 2.0);
    transient_coupling.thermal.initial_temperature.assign(
        box.mesh.nodes.size(), 300.0);
    nova::physics::TransientElectrothermal3DOptions transient_options;
    transient_options.thermal_time.time_step = 0.02;
    transient_options.thermal_time.steps = 2U;
    const auto transient_solution =
        nova::physics::solve_transient_electrothermal_3d(
            transient_coupling, transient_options);
    if (!transient_solution.converged ||
        transient_solution.accepted_steps.size() != 2U ||
        std::abs(transient_solution.last_electrical.total_joule_power_watts - 1.0) > 1e-8 ||
        std::abs(transient_solution.accepted_steps.back().heat_budget.balance_residual)
             > 1e-7)
        return 27;

    // The adaptive time-step controller must link from the installed SDK,
    // preserve discrete energy accounting and expose accepted-step telemetry.
    nova::physics::AdaptiveTransientElectrothermal3DOptions adaptive_thermal;
    adaptive_thermal.coupling = transient_options;
    adaptive_thermal.end_time = 0.02;
    adaptive_thermal.initial_time_step = 0.02;
    adaptive_thermal.maximum_time_step = 0.02;
    adaptive_thermal.absolute_error_kelvin = 0.1;
    const auto adaptive_solution =
        nova::physics::solve_adaptive_transient_electrothermal_3d(
            transient_coupling, adaptive_thermal);
    if (!adaptive_solution.completed ||
        adaptive_solution.accepted_steps.size() != 1U ||
        adaptive_solution.rejected_steps != 0U ||
        std::abs(adaptive_solution.final_time - 0.02) > 1e-12 ||
        std::abs(adaptive_solution.accepted_steps.front()
                     .heat_budget.balance_residual) > 1e-7)
        return 28;

    // Exercise the installed Neumann-boundary residual overload and link symbol.
    std::vector<NeumannResidualFace3D> natural_faces;
    for (const auto& face : box.boundary_faces) {
        const bool on_y0 = problem.mesh.nodes[face[0]].y == 0.0 &&
            problem.mesh.nodes[face[1]].y == 0.0 &&
            problem.mesh.nodes[face[2]].y == 0.0;
        if (on_y0) {
            natural_faces.push_back({
                face, [](double, double, double) { return 0.0; }
            });
        }
    }
    const std::vector<double> unit_conductivity(box.mesh.elements.size(), 1.0);
    const auto residual = estimate_diffusion_residual_error_3d(
        box.mesh, three_d.values, unit_conductivity,
        [](double, double, double) { return 0.0; }, natural_faces);
    if (natural_faces.empty() ||
        residual.boundary_residual_faces.size() != natural_faces.size() ||
        residual.boundary_residual_errors.size() != box.mesh.elements.size() ||
        residual.total_error > 1.0e-8) return 14;

    // Exercise installed adaptive 3D diffusion as a public linked SDK symbol.
    // A nonzero manufactured volume source guarantees a nontrivial mark/refine.
    ScalarDiffusion3DProblem adaptive_problem = problem;
    adaptive_problem.diffusivity = {};
    adaptive_problem.cell_diffusivity = unit_conductivity;
    adaptive_problem.source = [](double, double, double) { return 2.0; };
    for (auto& condition : adaptive_problem.dirichlet_faces) {
        condition.value = [](double, double, double) { return 0.0; };
    }
    const auto adaptive =
        solve_adaptive_diffusion_step_3d(adaptive_problem, 0.5);
    if (!adaptive.coarse_solution.linear_solver.converged ||
        !adaptive.refined_solution.linear_solver.converged ||
        adaptive.adaptation.marked_elements.empty() ||
        adaptive.refined_problem.mesh.nodes.size() <=
            adaptive_problem.mesh.nodes.size() ||
        adaptive.transferred_p1_values.size() !=
            adaptive.refined_problem.mesh.nodes.size() ||
        adaptive.refined_indicators.element_errors.size() !=
            adaptive.refined_problem.mesh.elements.size() ||
        !std::isfinite(adaptive.refined_indicators.total_error)) return 20;

    // Exercise the installed experimental electrostatics API and exported symbol.
    nova::physics::Electrostatics3DProblem electrostatic;
    electrostatic.mesh = box.mesh;
    electrostatic.permittivity = [](double, double, double) { return 2.0; };
    std::vector<nova::fem::BoundaryFace3D> electrode_faces;
    std::vector<int> electrode_ids;
    for (const auto& face : box.boundary_faces) {
        const bool left = electrostatic.mesh.nodes[face[0]].x == 0.0 &&
            electrostatic.mesh.nodes[face[1]].x == 0.0 &&
            electrostatic.mesh.nodes[face[2]].x == 0.0;
        const bool right = electrostatic.mesh.nodes[face[0]].x == 1.0 &&
            electrostatic.mesh.nodes[face[1]].x == 1.0 &&
            electrostatic.mesh.nodes[face[2]].x == 1.0;
        if (left || right) {
            electrostatic.potential_faces.push_back({face,
                [](double px, double, double) { return 10.0 * px; }});
            electrode_faces.push_back(face);
            electrode_ids.push_back(right ? 2 : 1);
        }
    }
    const auto capacitor = nova::physics::solve_electrostatics_3d(electrostatic);
    if (!capacitor.linear_solver.converged ||
        capacitor.potential.size() != electrostatic.mesh.nodes.size() ||
        std::abs(capacitor.stored_energy_joules - 100.0) > 1.0e-8) return 5;
    for (std::size_t i = 0; i < capacitor.potential.size(); ++i) {
        if (std::abs(capacitor.potential[i] - 10.0 * electrostatic.mesh.nodes[i].x) > 1.0e-9)
            return 6;
    }
    // Verify the installed fixed-electrode charge API and exported symbol.
    const auto electrode_charges = nova::physics::evaluate_electrode_charges_3d(
        electrostatic, capacitor, electrode_faces, electrode_ids);
    if (electrode_charges.size() != 2 ||
        electrode_charges[0].electrode_id != 1 ||
        electrode_charges[1].electrode_id != 2 ||
        std::abs(electrode_charges[0].conductor_charge_coulombs + 20.0) > 1.0e-8 ||
        std::abs(electrode_charges[1].conductor_charge_coulombs - 20.0) > 1.0e-8 ||
        std::abs(electrode_charges[1].area_square_metres - 1.0) > 1.0e-12) return 10;
    // The installed SDK also exposes the discrete P1 reaction charge,
    // independent of raw boundary-displacement surface integration.
    const auto reactions = nova::physics::evaluate_discrete_electrode_reactions_3d(
        electrostatic, capacitor, electrode_faces, electrode_ids);
    if (reactions.size() != 2 ||
        std::abs(reactions[0].reaction_charge_coulombs + 20.0) > 1.0e-8 ||
        std::abs(reactions[1].reaction_charge_coulombs - 20.0) > 1.0e-8)
        return 17;
    const double installed_capacitance =
        electrode_charges[1].conductor_charge_coulombs /
        (electrode_charges[1].potential_volts - electrode_charges[0].potential_volts);
    if (std::abs(installed_capacitance - 2.0) > 1.0e-9) return 11;
    // Verify the newly installed generic Gmsh importer as an SDK consumer.
    // This minimal MSH 4.1 file has sparse-compatible block records and no
    // physical groups; topology and exterior faces are reconstructed by NOVA.
    const std::string msh41 = R"msh($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1 2 3 4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
$EndElements
)msh";
    std::istringstream mesh_stream(msh41);
    const auto imported = nova::io::read_gmsh_tetrahedra(mesh_stream);
    if (imported.mesh.nodes.size() != 4 ||
        imported.mesh.elements.size() != 1 ||
        imported.volume_physical_tags.size() != 1 ||
        imported.volume_physical_tags[0] != 0 ||
        imported.exterior_faces.size() != 4) return 7;

    // Verify the installed multi-terminal capacitance-matrix API and symbol.
    std::vector<nova::physics::ElectrodeBoundaryGroup3D> electrode_groups{
        {1, {}}, {2, {}}
    };
    for (const auto& face : box.boundary_faces) {
        const bool left = electrostatic.mesh.nodes[face[0]].x == 0.0 &&
            electrostatic.mesh.nodes[face[1]].x == 0.0 &&
            electrostatic.mesh.nodes[face[2]].x == 0.0;
        const bool right = electrostatic.mesh.nodes[face[0]].x == 1.0 &&
            electrostatic.mesh.nodes[face[1]].x == 1.0 &&
            electrostatic.mesh.nodes[face[2]].x == 1.0;
        if (left) electrode_groups[0].faces.push_back(face);
        if (right) electrode_groups[1].faces.push_back(face);
    }
    const auto capacitance_matrix =
        nova::physics::extract_capacitance_matrix_3d(electrostatic, electrode_groups);
    if (capacitance_matrix.electrode_ids != std::vector<int>{1, 2} ||
        capacitance_matrix.matrix_farads.size() != 4 ||
        std::abs(capacitance_matrix.matrix_farads[0] - 2.0) > 1.0e-8 ||
        std::abs(capacitance_matrix.matrix_farads[1] + 2.0) > 1.0e-8 ||
        std::abs(capacitance_matrix.matrix_farads[2] + 2.0) > 1.0e-8 ||
        std::abs(capacitance_matrix.matrix_farads[3] - 2.0) > 1.0e-8 ||
        capacitance_matrix.max_reciprocity_defect_farads > 1.0e-8) return 13;

    // Verify the installed FEM-reaction capacitance extractor as an
    // independent public symbol, preserving the raw flux API.
    const auto discrete_capacitance =
        nova::physics::extract_discrete_capacitance_matrix_3d(
            electrostatic, electrode_groups);
    if (discrete_capacitance.electrode_ids !=
            std::vector<int>({1, 2}) ||
        discrete_capacitance.matrix_farads.size() != 4U ||
        std::abs(discrete_capacitance.matrix_farads[0] - 2.0) > 1.0e-8 ||
        std::abs(discrete_capacitance.matrix_farads[1] + 2.0) > 1.0e-8 ||
        discrete_capacitance.max_reciprocity_defect_farads > 1.0e-8)
        return 18;

    // Check installed floating-electrode terminal algebra and exported symbol.
    const auto floating_voltages =
        nova::physics::solve_floating_electrode_potentials_3d(
            capacitance_matrix,
            std::vector<nova::physics::FixedElectrodePotential3D>{{1, 0.0}},
            std::vector<nova::physics::FloatingElectrodeCharge3D>{{2, 20.0}});
    if (floating_voltages.potential_volts.size() != 2 ||
        std::abs(floating_voltages.potential_volts[1] - 10.0) > 1.0e-8 ||
        std::abs(floating_voltages.estimated_charge_coulombs[0] + 20.0) > 1.0e-8 ||
        floating_voltages.max_floating_charge_residual_coulombs > 1.0e-8)
        return 15;

    // Installed SDK symbol: terminal inference followed by a new FEM field solve.
    const auto reconstructed =
        nova::physics::reconstruct_floating_electrostatic_field_3d(
            electrostatic, electrode_groups,
            std::vector<nova::physics::FixedElectrodePotential3D>{{1, 0.0}},
            std::vector<nova::physics::FloatingElectrodeCharge3D>{{2, 20.0}});
    if (!reconstructed.field.linear_solver.converged ||
        reconstructed.field.potential.size() != box.mesh.nodes.size() ||
        reconstructed.electrode_charges.size() != 2 ||
        reconstructed.max_floating_charge_mismatch_coulombs > 1.0e-8 ||
        std::abs(reconstructed.field.stored_energy_joules - 100.0) > 1.0e-7)
        return 16;

    // Installed discrete Schur system: Q_high=20 C -> V_high=10 V
    // -> one reconstructed FEM field and consistent Galerkin electrode Q.
    const auto schur = nova::physics::solve_floating_conductor_schur_3d(
        electrostatic, electrode_groups,
        std::vector<nova::physics::FixedElectrodePotential3D>{{1, 0.0}},
        std::vector<nova::physics::FloatingElectrodeCharge3D>{{2, 20.0}});
    if (!schur.field.linear_solver.converged ||
        schur.discrete_reactions.size() != 2 ||
        schur.boundary_flux_charges.size() != 2 ||
        std::abs(schur.terminals.potential_volts[1] - 10.0) > 1.0e-8 ||
        std::abs(schur.discrete_reactions[1].reaction_charge_coulombs -
                 20.0) > 1.0e-8 ||
        schur.max_floating_reaction_mismatch_coulombs > 1.0e-8 ||
        std::abs(schur.field.stored_energy_joules - 100.0) > 1.0e-7)
        return 19;

    // Verify the installed SDK exports the versioned thermal checkpoint API.
    const std::vector<double> saved_temperature(box.mesh.nodes.size(), 301.25);
    std::stringstream checkpoint_stream;
    nova::io::write_transient_heat_checkpoint_3d(
        checkpoint_stream, box.mesh, 1.75, saved_temperature);
    std::istringstream checkpoint_input(checkpoint_stream.str());
    const auto restored = nova::io::read_transient_heat_checkpoint_3d(
        checkpoint_input, box.mesh);
    if (restored.time != 1.75 || restored.temperature != saved_temperature) return 12;
    return 0;
}
