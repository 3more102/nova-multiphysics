#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/energy_balance_3d.hpp"

#include <vector>

namespace nova::physics {

// Experimental stationary ONE-WAY DC electrical -> thermal FEM coupling.
// -div(sigma grad(V)) = impressed current density [A/m^3].
// -div(k grad(T)) = q_base + grad(V) dot sigma grad(V).
// Both physics use precisely the same tetrahedral mesh and ordering.
// There is no temperature-dependent electrical conductivity in this API.
struct Electrothermal3DProblem {
    fem::ScalarDiffusion3DProblem electrical; // sigma [S/m], V [V].
    fem::ScalarDiffusion3DProblem thermal;    // k [W/(m K)], T [K].
};

struct Electrothermal3DResult {
    fem::ScalarDiffusion3DResult electrical_solution;
    std::vector<fem::Point3D> electric_field;  // -grad(V), [V/m].
    std::vector<fem::Point3D> current_density; // -sigma grad(V), [A/m^2].
    std::vector<double> joule_density;         // Cell-average E dot J [W/m^3].
    double total_joule_power_watts{};
    fem::ScalarDiffusion3DProblem heated_thermal_problem;
    fem::ScalarDiffusion3DResult thermal_solution;
    fem::DiffusionEnergyBalance3D thermal_balance;
};

[[nodiscard]] Electrothermal3DResult solve_one_way_electrothermal_3d(
    const Electrothermal3DProblem& problem,
    const linalg::SolverOptions& electrical_options = {},
    const linalg::SolverOptions& thermal_options = {}
);

} // namespace nova::physics
