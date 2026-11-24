#include "nova/physics/photothermal_3d.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace nova::physics {
namespace {
double checked_double(long double value) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max()))
        throw std::overflow_error("Photothermal 3D heat source exceeds double range.");
    return static_cast<double>(value);
}
}

Photothermal3DResult solve_one_way_photothermal_3d(
    const Photothermal3DProblem& problem,
    const linalg::GmresOptions& optical_options,
    const linalg::SolverOptions& thermal_options) {
    const auto& mesh = problem.thermal.mesh;
    const auto n = mesh.elements.size();
    const auto& base = problem.thermal.cell_source;
    if (problem.optical_materials.size() != n ||
        problem.optical_cell_sources.size() != n ||
        (!base.empty() && base.size() != n))
        throw std::invalid_argument("Photothermal 3D cell counts must match the mesh.");
    for (double q : base)
        if (!std::isfinite(q))
            throw std::invalid_argument("Photothermal 3D base source must be finite.");

    const auto system = assemble_maxwell_3d(mesh, problem.optical_materials,
        problem.vacuum_wavenumber, problem.pec_on_exterior);
    const auto load = assemble_maxwell_cell_load_3d(
        mesh, problem.optical_cell_sources, system);
    Photothermal3DResult result;
    result.optical_solution = solve_maxwell_3d(system, load, optical_options);
    if (!result.optical_solution.diagnostics.converged)
        throw std::runtime_error("Photothermal 3D Maxwell solve did not converge.");
    result.optical_edges = system.edges;
    result.absorbed_power_per_cell = integrate_maxwell_absorption_3d(
        mesh, problem.optical_materials, problem.vacuum_wavenumber,
        system, result.optical_solution.edge_field);
    result.heated_thermal_problem = problem.thermal;
    auto& sources = result.heated_thermal_problem.cell_source;
    if (sources.empty()) sources.assign(n, 0.0);
    result.absorbed_power_density.reserve(n);
    long double total = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        const double power = result.absorbed_power_per_cell[i];
        const double volume = fem::make_nedelec_tetrahedron_3d(mesh, i).volume;
        if (!std::isfinite(power) || power < 0.0 ||
            !std::isfinite(volume) || volume <= 0.0)
            throw std::runtime_error("Photothermal 3D invalid absorption or volume.");
        const auto density = static_cast<long double>(power) / volume;
        result.absorbed_power_density.push_back(checked_double(density));
        sources[i] = checked_double(static_cast<long double>(sources[i]) + density);
        total += power;
    }
    result.total_absorbed_power_watts = checked_double(total);
    result.thermal_solution = fem::solve_scalar_diffusion_3d(
        result.heated_thermal_problem, thermal_options);
    if (!result.thermal_solution.linear_solver.converged)
        throw std::runtime_error("Photothermal 3D thermal solve did not converge.");
    result.thermal_balance = fem::evaluate_diffusion_energy_balance_3d(
        result.heated_thermal_problem, result.thermal_solution.values);
    return result;
}

} // namespace nova::physics
