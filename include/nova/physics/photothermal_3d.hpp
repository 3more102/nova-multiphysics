#pragma once

#include "nova/physics/maxwell_3d.hpp"
#include "nova/fem/energy_balance_3d.hpp"

namespace nova::physics {

// Stationary one-way optical absorption -> heat conduction on a single mesh.
// Thermal mesh defines the geometry and ordering for both solves.
// Optical coefficients are frozen: no thermo-optic feedback or nonlinear optics.
struct Photothermal3DProblem {
    fem::ScalarDiffusion3DProblem thermal;
    std::vector<MaxwellCellMaterial3D> optical_materials;
    std::vector<ComplexVector3D> optical_cell_sources;
    double vacuum_wavenumber{};
    bool pec_on_exterior{true};
};

struct Photothermal3DResult {
    MaxwellSolution3D optical_solution;
    std::vector<fem::Edge3D> optical_edges;
    std::vector<double> absorbed_power_per_cell; // W
    std::vector<double> absorbed_power_density; // W/m^3
    double total_absorbed_power_watts{};
    fem::ScalarDiffusion3DProblem heated_thermal_problem;
    fem::ScalarDiffusion3DResult thermal_solution;
    fem::DiffusionEnergyBalance3D thermal_balance;
};

// Uses exact integrated dielectric absorption, divided by each cell volume,
// to preserve deposited power in the P1 thermal weak source. Existing base
// sources are additive. Throws if either linear solve fails to converge.
[[nodiscard]] Photothermal3DResult solve_one_way_photothermal_3d(
    const Photothermal3DProblem& problem,
    const linalg::GmresOptions& optical_options = {},
    const linalg::SolverOptions& thermal_options = {});

} // namespace nova::physics
