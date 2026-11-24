#pragma once

#include "nova/physics/photothermal_3d.hpp"

#include <complex>
#include <cstddef>
#include <vector>

namespace nova::physics {

// Cellwise linear thermo-optic law:
// epsilon_r(T) = epsilon_r,reference + d_epsilon_r_dT * (T_cell - T_reference).
// T_cell is the arithmetic mean of the four P1 vertex temperatures.
struct ThermoOptic3DProblem {
    Photothermal3DProblem photothermal;
    std::vector<std::complex<double>> d_relative_permittivity_d_temperature;
    double reference_temperature_kelvin{300.0};
    double relaxation_factor{0.5};
    double temperature_tolerance_kelvin{1e-6};
    double relative_absorbed_power_tolerance{1e-6};
    std::size_t max_iterations{50};
};

struct ThermoOptic3DIteration {
    std::size_t iteration{};
    double max_temperature_residual_kelvin{};
    double relative_absorbed_power_change{};
    double total_absorbed_power_watts{};
};

struct ThermoOptic3DResult {
    Photothermal3DResult photothermal;
    // Permittivity used by the final Maxwell solve, in mesh cell order.
    std::vector<MaxwellCellMaterial3D> optical_materials;
    std::vector<ThermoOptic3DIteration> history;
    bool converged{};
};

// Iterates frequency-domain Maxwell absorption and stationary P1 heat
// conduction to a self-consistent temperature-dependent permittivity.
// Returns converged=false with the last finite iterate if max_iterations is
// reached. Invalid inputs or a failed inner linear solve throw.
[[nodiscard]] ThermoOptic3DResult solve_thermo_optic_3d(
    const ThermoOptic3DProblem& problem,
    const linalg::GmresOptions& optical_options = {},
    const linalg::SolverOptions& thermal_options = {});

} // namespace nova::physics
