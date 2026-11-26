#include "nova/physics/thermo_optic_3d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace nova::physics {
namespace {

bool finite(std::complex<double> value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void validate_problem(const ThermoOptic3DProblem& problem) {
    const auto& coupled = problem.photothermal;
    const auto& mesh = coupled.thermal.mesh;
    mesh.validate();
    const auto cells = mesh.elements.size();
    if (cells == 0 ||
        coupled.optical_materials.size() != cells ||
        coupled.optical_cell_sources.size() != cells ||
        problem.d_relative_permittivity_d_temperature.size() != cells ||
        (!coupled.thermal.cell_source.empty() &&
         coupled.thermal.cell_source.size() != cells)) {
        throw std::invalid_argument("Thermo-optic 3D cell data must match the tetrahedral mesh.");
    }
    if (!std::isfinite(coupled.vacuum_wavenumber) ||
        coupled.vacuum_wavenumber <= 0.0 ||
        !std::isfinite(problem.reference_temperature_kelvin) ||
        problem.reference_temperature_kelvin <= 0.0 ||
        !std::isfinite(problem.relaxation_factor) ||
        problem.relaxation_factor <= 0.0 || problem.relaxation_factor > 1.0 ||
        !std::isfinite(problem.temperature_tolerance_kelvin) ||
        problem.temperature_tolerance_kelvin <= 0.0 ||
        !std::isfinite(problem.relative_absorbed_power_tolerance) ||
        problem.relative_absorbed_power_tolerance <= 0.0 ||
        !std::isfinite(problem.permittivity_tolerance) ||
        problem.permittivity_tolerance <= 0.0 ||
        problem.max_iterations == 0) {
        throw std::invalid_argument("Invalid thermo-optic 3D iteration controls.");
    }
    for (std::size_t i = 0; i < cells; ++i) {
        const auto epsilon = coupled.optical_materials[i].relative_permittivity;
        if (!finite(epsilon) || epsilon.imag() < 0.0 ||
            !finite(problem.d_relative_permittivity_d_temperature[i])) {
            throw std::invalid_argument(
                "Thermo-optic 3D materials must have finite passive reference permittivity and finite slope.");
        }
    }
}

} // namespace

ThermoOptic3DResult solve_thermo_optic_3d(
    const ThermoOptic3DProblem& problem,
    const linalg::GmresOptions& optical_options,
    const linalg::SolverOptions& thermal_options) {
    validate_problem(problem);

    const auto& base = problem.photothermal;
    const auto& mesh = base.thermal.mesh;

    // Start from the stationary thermal field without optical absorption. This
    // honours the user's base source and boundary conditions from iteration 0.
    auto baseline = fem::solve_scalar_diffusion_3d(base.thermal, thermal_options);
    if (!baseline.linear_solver.converged) {
        throw std::runtime_error("Thermo-optic 3D baseline thermal solve did not converge.");
    }
    auto current_temperature = std::move(baseline.values);
    if (current_temperature.size() != mesh.nodes.size()) {
        throw std::runtime_error("Thermo-optic 3D baseline field size mismatch.");
    }

    ThermoOptic3DResult result;
    double previous_absorbed_power = 0.0;

    for (std::size_t iteration = 0; iteration < problem.max_iterations; ++iteration) {
        auto materials = base.optical_materials;
        for (std::size_t cell = 0; cell < mesh.elements.size(); ++cell) {
            double cell_temperature = 0.0;
            for (const auto node : mesh.elements[cell]) {
                cell_temperature += 0.25 * current_temperature[node];
            }
            const double delta_temperature =
                cell_temperature - problem.reference_temperature_kelvin;
            const auto epsilon = base.optical_materials[cell].relative_permittivity +
                problem.d_relative_permittivity_d_temperature[cell] * delta_temperature;
            if (!finite(epsilon) || epsilon.imag() < 0.0) {
                throw std::invalid_argument(
                    "Thermo-optic 3D update produced nonfinite or active permittivity.");
            }
            materials[cell].relative_permittivity = epsilon;
        }

        auto iteration_problem = base;
        iteration_problem.optical_materials = materials;
        auto coupled = solve_one_way_photothermal_3d(
            iteration_problem, optical_options, thermal_options);

        const auto& next_temperature = coupled.thermal_solution.values;
        if (next_temperature.size() != current_temperature.size()) {
            throw std::runtime_error("Thermo-optic 3D thermal field size changed between iterations.");
        }

        double max_temperature_residual = 0.0;
        for (std::size_t node = 0; node < current_temperature.size(); ++node) {
            if (!std::isfinite(next_temperature[node])) {
                throw std::runtime_error("Thermo-optic 3D produced a nonfinite temperature.");
            }
            max_temperature_residual = std::max(
                max_temperature_residual,
                std::abs(next_temperature[node] - current_temperature[node]));
        }

        // The returned thermal field is not the input iterate used to
        // assemble Maxwell. Check constitutive closure explicitly: a small
        // temperature error alone need not imply a small permittivity error.
        double max_permittivity_residual = 0.0;
        bool output_materials_passive = true;
        for (std::size_t cell = 0; cell < mesh.elements.size(); ++cell) {
            double output_cell_temperature = 0.0;
            for (const auto node : mesh.elements[cell]) {
                output_cell_temperature += 0.25 * next_temperature[node];
            }
            const auto output_epsilon = base.optical_materials[cell].relative_permittivity +
                problem.d_relative_permittivity_d_temperature[cell] *
                (output_cell_temperature - problem.reference_temperature_kelvin);
            const double defect = std::abs(
                output_epsilon - materials[cell].relative_permittivity);
            if (!finite(output_epsilon) || !std::isfinite(defect)) {
                throw std::runtime_error("Thermo-optic 3D material closure is nonfinite.");
            }
            output_materials_passive = output_materials_passive && output_epsilon.imag() >= 0.0;
            max_permittivity_residual = std::max(max_permittivity_residual, defect);
        }

        const double absorbed_power = coupled.total_absorbed_power_watts;
        const double power_scale = std::max({
            std::abs(absorbed_power), std::abs(previous_absorbed_power), 1e-30});
        const double relative_power_change =
            std::abs(absorbed_power - previous_absorbed_power) / power_scale;

        result.history.push_back({
            iteration + 1,
            max_temperature_residual,
            relative_power_change,
            absorbed_power,
            max_permittivity_residual
        });

        if (max_temperature_residual <= problem.temperature_tolerance_kelvin &&
            relative_power_change <= problem.relative_absorbed_power_tolerance &&
            max_permittivity_residual <= problem.permittivity_tolerance &&
            output_materials_passive) {
            result.photothermal = std::move(coupled);
            result.optical_materials = std::move(materials);
            result.converged = true;
            return result;
        }

        // Consume the next field before moving the coupled result into the
        // persistent return object; moving first would invalidate this reference.
        for (std::size_t node = 0; node < current_temperature.size(); ++node) {
            current_temperature[node] += problem.relaxation_factor *
                (next_temperature[node] - current_temperature[node]);
            if (!std::isfinite(current_temperature[node])) {
                throw std::runtime_error("Thermo-optic 3D relaxation produced a nonfinite temperature.");
            }
        }
        previous_absorbed_power = absorbed_power;
        result.photothermal = std::move(coupled);
        result.optical_materials = std::move(materials);
    }

    return result;
}

} // namespace nova::physics
