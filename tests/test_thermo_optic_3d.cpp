#include "test_framework.hpp"
#include "nova/physics/thermo_optic_3d.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace {
using namespace nova::fem;
using namespace nova::physics;

Photothermal3DProblem plate(double loss, double amplitude) {
    const auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    Photothermal3DProblem problem;
    problem.thermal.mesh = box.mesh;
    problem.thermal.cell_diffusivity.assign(box.mesh.elements.size(), 4.0);
    problem.thermal.cell_source.assign(box.mesh.elements.size(), 7.0);
    problem.thermal.source = [](double,double,double) { return 3.0; };
    problem.vacuum_wavenumber = 2.3;
    problem.pec_on_exterior = false;

    const std::complex<double> epsilon(2.0, loss);
    problem.optical_materials.assign(box.mesh.elements.size(), {epsilon, 1.0});
    ComplexVector3D source{};
    source[0] = -problem.vacuum_wavenumber * problem.vacuum_wavenumber *
        epsilon * amplitude;
    problem.optical_cell_sources.assign(box.mesh.elements.size(), source);

    for (const auto& face : box.boundary_faces) {
        const auto on_x = [&](double x) {
            return std::all_of(face.begin(), face.end(),
                [&](std::size_t node) { return box.mesh.nodes[node].x == x; });
        };
        if (on_x(0.0) || on_x(1.0)) {
            problem.thermal.dirichlet_faces.push_back(
                {face, [](double,double,double) { return 300.0; }});
        }
    }
    return problem;
}

ThermoOptic3DProblem thermo_optic_problem(const Photothermal3DProblem& base) {
    ThermoOptic3DProblem problem;
    problem.photothermal = base;
    problem.d_relative_permittivity_d_temperature.assign(
        base.thermal.mesh.elements.size(), {0.0, 0.0});
    problem.reference_temperature_kelvin = 300.0;
    problem.relaxation_factor = 1.0;
    problem.temperature_tolerance_kelvin = 1e-8;
    problem.relative_absorbed_power_tolerance = 1e-8;
    problem.max_iterations = 50;
    return problem;
}

nova::linalg::GmresOptions accurate_optical_options() {
    nova::linalg::GmresOptions options;
    options.restart = 100;
    options.solver.max_iterations = 1000;
    options.solver.relative_tolerance = 1e-12;
    return options;
}
}

NOVA_TEST(thermo_optic_zero_slope_matches_one_way_photothermal) {
    const auto base = plate(0.6, 100.0);
    const auto reference = solve_one_way_photothermal_3d(base, accurate_optical_options());
    const auto problem = thermo_optic_problem(base);
    const auto result = solve_thermo_optic_3d(problem, accurate_optical_options());

    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.history.size() == 2);
    NOVA_REQUIRE(result.photothermal.optical_solution.diagnostics.converged);
    NOVA_REQUIRE(result.photothermal.thermal_solution.linear_solver.converged);
    NOVA_REQUIRE_NEAR(result.photothermal.total_absorbed_power_watts,
                      reference.total_absorbed_power_watts, 1e-9);
    NOVA_REQUIRE_NEAR(result.photothermal.thermal_balance.signed_balance_defect, 0.0, 1e-6);
    NOVA_REQUIRE(result.photothermal.thermal_solution.values.size() ==
                 reference.thermal_solution.values.size());
    for (std::size_t i = 0; i < reference.thermal_solution.values.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.photothermal.thermal_solution.values[i],
                          reference.thermal_solution.values[i], 1e-9);
    }
    for (std::size_t i = 0; i < result.optical_materials.size(); ++i) {
        NOVA_REQUIRE_NEAR(result.optical_materials[i].relative_permittivity.real(), 2.0, 1e-14);
        NOVA_REQUIRE_NEAR(result.optical_materials[i].relative_permittivity.imag(), 0.6, 1e-14);
    }
}

NOVA_TEST(thermo_optic_temperature_feedback_converges_with_energy_balance) {
    const auto base = plate(0.6, 100.0);
    auto problem = thermo_optic_problem(base);
    problem.d_relative_permittivity_d_temperature.assign(
        base.thermal.mesh.elements.size(), {1e-3, 2e-4});
    problem.relaxation_factor = 0.5;
    problem.temperature_tolerance_kelvin = 1e-7;
    problem.relative_absorbed_power_tolerance = 1e-7;
    problem.max_iterations = 80;

    const auto result = solve_thermo_optic_3d(problem, accurate_optical_options());
    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.history.size() > 2);
    NOVA_REQUIRE(result.history.back().max_temperature_residual_kelvin <=
                 problem.temperature_tolerance_kelvin);
    NOVA_REQUIRE(result.history.back().relative_absorbed_power_change <=
                 problem.relative_absorbed_power_tolerance);
    NOVA_REQUIRE(result.photothermal.total_absorbed_power_watts > 0.0);
    NOVA_REQUIRE_NEAR(result.photothermal.thermal_balance.signed_balance_defect, 0.0, 1e-6);

    const auto& temperatures = result.photothermal.thermal_solution.values;
    for (std::size_t cell = 0; cell < base.thermal.mesh.elements.size(); ++cell) {
        double mean_temperature = 0.0;
        for (const auto node : base.thermal.mesh.elements[cell]) {
            mean_temperature += 0.25 * temperatures[node];
        }
        const auto expected = base.optical_materials[cell].relative_permittivity +
            problem.d_relative_permittivity_d_temperature[cell] *
                (mean_temperature - problem.reference_temperature_kelvin);
        NOVA_REQUIRE_NEAR(result.optical_materials[cell].relative_permittivity.real(),
                          expected.real(), 1e-9);
        NOVA_REQUIRE_NEAR(result.optical_materials[cell].relative_permittivity.imag(),
                          expected.imag(), 1e-9);
    }
}

NOVA_TEST(thermo_optic_reports_nonconvergence_and_rejects_bad_cell_data) {
    const auto base = plate(0.6, 100.0);
    auto problem = thermo_optic_problem(base);
    problem.max_iterations = 1;
    const auto incomplete = solve_thermo_optic_3d(problem, accurate_optical_options());
    NOVA_REQUIRE(!incomplete.converged);
    NOVA_REQUIRE(incomplete.history.size() == 1);
    NOVA_REQUIRE(incomplete.history.front().relative_absorbed_power_change > 0.0);

    problem.d_relative_permittivity_d_temperature.pop_back();
    bool rejected = false;
    try {
        (void)solve_thermo_optic_3d(problem, accurate_optical_options());
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    NOVA_REQUIRE(rejected);
}


NOVA_TEST(thermo_optic_manufactured_nonzero_temperature_shift_matches_exact_absorption) {
    const auto box = make_box_mesh_3d(0,1,0,1,0,1,2,2,2);
    Photothermal3DProblem base;
    base.thermal.mesh = box.mesh;
    base.thermal.cell_diffusivity.assign(box.mesh.elements.size(), 4.0);
    for (std::size_t node = 0; node < box.mesh.nodes.size(); ++node) {
        base.thermal.dirichlet_constraints.push_back({node, 320.0});
    }

    constexpr double k0 = 2.3;
    constexpr double amplitude = 100.0;
    constexpr double c = 299792458.0;
    constexpr double eps0 = 8.8541878128e-12;
    const std::complex<double> epsilon_reference(2.0, 0.6);
    const std::complex<double> slope(1e-3, 2e-4);
    const auto epsilon_target = epsilon_reference + slope * 20.0;

    base.vacuum_wavenumber = k0;
    base.pec_on_exterior = false;
    base.optical_materials.assign(box.mesh.elements.size(), {epsilon_reference, 1.0});
    ComplexVector3D source{};
    // Manufacture a uniform E=(100,0,0) field for the updated permittivity:
    // curl(curl(E)) - k0^2 * epsilon_target * E = source.
    source[0] = -k0 * k0 * epsilon_target * amplitude;
    base.optical_cell_sources.assign(box.mesh.elements.size(), source);

    ThermoOptic3DProblem problem;
    problem.photothermal = base;
    problem.d_relative_permittivity_d_temperature.assign(
        box.mesh.elements.size(), slope);
    problem.reference_temperature_kelvin = 300.0;
    problem.relaxation_factor = 1.0;
    problem.temperature_tolerance_kelvin = 1e-10;
    problem.relative_absorbed_power_tolerance = 1e-10;
    problem.max_iterations = 4;

    const auto result = solve_thermo_optic_3d(problem, accurate_optical_options());
    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.history.size() == 2);
    NOVA_REQUIRE(result.history.front().relative_absorbed_power_change > 0.99);
    NOVA_REQUIRE_NEAR(result.history.back().relative_absorbed_power_change, 0.0, 1e-12);

    const double expected_density =
        0.5 * k0 * c * eps0 * epsilon_target.imag() * amplitude * amplitude;
    NOVA_REQUIRE_NEAR(result.photothermal.total_absorbed_power_watts,
                      expected_density, 1e-7);
    NOVA_REQUIRE_NEAR(result.photothermal.thermal_balance.generated_heat,
                      expected_density, 1e-7);
    NOVA_REQUIRE_NEAR(result.photothermal.thermal_balance.constraint_outward_heat,
                      expected_density, 1e-6);
    NOVA_REQUIRE_NEAR(result.photothermal.thermal_balance.signed_balance_defect,
                      0.0, 1e-6);

    NOVA_REQUIRE(result.photothermal.thermal_solution.values.size() == box.mesh.nodes.size());
    for (double temperature : result.photothermal.thermal_solution.values) {
        NOVA_REQUIRE_NEAR(temperature, 320.0, 1e-12);
    }
    for (std::size_t cell = 0; cell < box.mesh.elements.size(); ++cell) {
        const double volume = make_nedelec_tetrahedron_3d(box.mesh, cell).volume;
        NOVA_REQUIRE_NEAR(result.optical_materials[cell].relative_permittivity.real(),
                          epsilon_target.real(), 1e-12);
        NOVA_REQUIRE_NEAR(result.optical_materials[cell].relative_permittivity.imag(),
                          epsilon_target.imag(), 1e-12);
        NOVA_REQUIRE_NEAR(result.photothermal.absorbed_power_density[cell],
                          expected_density, 1e-7);
        NOVA_REQUIRE_NEAR(result.photothermal.absorbed_power_per_cell[cell],
                          expected_density * volume, 1e-7);
    }
}

// The optical field remains of ordinary magnitude; only dielectric loss
// is tiny. This isolates the outer convergence criterion from Krylov scaling.
NOVA_TEST(thermo_optic_weak_absorption_requires_relative_power_convergence) {
    const auto base = plate(1e-45, 100.0);
    auto problem = thermo_optic_problem(base);
    problem.max_iterations = 1;
    const auto incomplete = solve_thermo_optic_3d(problem, accurate_optical_options());

    const double power = incomplete.photothermal.total_absorbed_power_watts;
    NOVA_REQUIRE(power > 0.0);
    NOVA_REQUIRE(power < 1e-38);
    NOVA_REQUIRE(incomplete.history.size() == 1);
    NOVA_REQUIRE(incomplete.history.front().max_temperature_residual_kelvin <=
                 problem.temperature_tolerance_kelvin);
    NOVA_REQUIRE_NEAR(incomplete.history.front().relative_absorbed_power_change,
                      1.0, 1e-14);
    NOVA_REQUIRE(!incomplete.converged);

    problem.max_iterations = 2;
    const auto complete = solve_thermo_optic_3d(problem, accurate_optical_options());
    NOVA_REQUIRE(complete.converged);
    NOVA_REQUIRE(complete.history.size() == 2);
    NOVA_REQUIRE_NEAR(complete.history.back().relative_absorbed_power_change,
                      0.0, 1e-14);
    constexpr double c = 299792458.0;
    constexpr double eps0 = 8.8541878128e-12;
    const double expected = 0.5 * base.vacuum_wavenumber * c * eps0 *
        1e-45 * 100.0 * 100.0;
    NOVA_REQUIRE_NEAR(complete.photothermal.total_absorbed_power_watts / expected,
                      1.0, 1e-8);
}

NOVA_TEST(thermo_optic_exact_zero_absorption_has_zero_relative_change) {
    auto problem = thermo_optic_problem(plate(0.0, 100.0));
    problem.max_iterations = 1;
    const auto result = solve_thermo_optic_3d(problem, accurate_optical_options());
    NOVA_REQUIRE(result.converged);
    NOVA_REQUIRE(result.history.size() == 1);
    NOVA_REQUIRE(result.photothermal.total_absorbed_power_watts == 0.0);
    NOVA_REQUIRE(result.history.front().relative_absorbed_power_change == 0.0);
}
