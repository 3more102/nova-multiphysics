#pragma once

#include "nova/physics/electrode_charge_3d.hpp"

#include <span>
#include <vector>

namespace nova::physics {

// Known terminal voltages [V] and unknown-voltage terminal charges [C].
// Every electrode ID in the capacitance matrix must occur in exactly one
// of these two input collections.
struct FixedElectrodePotential3D {
    int electrode_id{};
    double potential_volts{};
};

struct FloatingElectrodeCharge3D {
    int electrode_id{};
    double charge_coulombs{};
};

struct FloatingElectrodePotentials3D {
    // Sorted in the same order as CapacitanceMatrix3D::electrode_ids.
    std::vector<int> electrode_ids;
    std::vector<double> potential_volts;
    std::vector<double> estimated_charge_coulombs;
    double max_floating_charge_residual_coulombs{};
};

// Solve C_FF V_F = Q_F - C_FK V_K for unknown terminal voltages.
// Uses the RAW, unsymmetrized estimated capacitance matrix. This performs
// terminal algebra only: it does not re-solve electrostatic FEM fields.
// A zero-reference Maxwell matrix can be singular if no terminal is fixed;
// a separately grounded external boundary may make it nonsingular.
// Requires a well-conditioned floating submatrix, with an approximately
// 1e-12 relative pivot threshold. Throws on invalid data or singularity.
// The extraction's charge-free and numerical-flux accuracy limits apply.
[[nodiscard]] FloatingElectrodePotentials3D solve_floating_electrode_potentials_3d(
    const CapacitanceMatrix3D& capacitance,
    std::span<const FixedElectrodePotential3D> fixed_electrodes,
    std::span<const FloatingElectrodeCharge3D> floating_electrodes
);

} // namespace nova::physics
