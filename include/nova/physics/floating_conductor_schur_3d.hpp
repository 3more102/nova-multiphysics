#pragma once

#include "nova/physics/floating_electrodes_3d.hpp"

#include <span>
#include <vector>

namespace nova::physics {

// Experimental 3D P1 floating-conductor solve for linear dielectrics.
// Unit-electrode Dirichlet FEM solves implicitly eliminate interior unknowns,
// producing a discrete Schur terminal matrix C = K_ee-K_ei*K_ii^{-1}*K_ie.
// With impressed volume free charge or outward Neumann displacement, one
// additional zero-terminal-voltage FEM solve obtains the induced weak reaction
// charge q0. Unknown terminal voltages then obey
//   C_FF V_F = Q_F - q0_F - C_FK V_K.
// A final loaded FEM solve reconstructs spatial E and D fields.
//
// External reference constraints/faces are grounded at 0 V, matching the
// capacitance basis solves; original nonzero external Dirichlet values are
// NOT preserved. The terminal matrix remains the SOURCE-FREE response;
// reported terminal charges are C*V+q0, not C*V alone. The raw boundary-flux
// estimate is not equilibrated and may disagree with loaded weak reactions.
// A gauge-singular floating block is rejected. This is an implicit discrete
// Schur workflow, not a monolithic saddle-point formulation.
struct SchurFloatingConductor3DResult {
    CapacitanceMatrix3D discrete_capacitance;
    // Sorted by discrete_capacitance.electrode_ids; charges from impressed
    // loads with ALL electrode voltages held at zero [C].
    // Exactly zero for source-free models (no offset solve).
    std::vector<double> induced_terminal_charges_coulombs;
    // Estimated charges include q0, and floating residuals refer to the
    // original requested physical charges rather than shifted RHS charges.
    FloatingElectrodePotentials3D terminals;
    Electrostatics3DResult field;
    std::vector<DiscreteElectrodeReaction3D> discrete_reactions;
    // Independent boundary-displacement estimate (not equilibrated).
    std::vector<ElectrodeCharge3D> boundary_flux_charges;
    // Max |Q_F^weak - Q_F^requested| over floating terminals [C].
    double max_floating_reaction_mismatch_coulombs{};
    // Max |Q_g^weak - Q_g^surface| over all terminals [C].
    // This is a mesh-sensitive cross-diagnostic, not a certified error.
    double max_flux_reaction_disagreement_coulombs{};
};

[[nodiscard]] SchurFloatingConductor3DResult solve_floating_conductor_schur_3d(
    const Electrostatics3DProblem& dielectric_problem,
    std::span<const ElectrodeBoundaryGroup3D> electrode_groups,
    std::span<const FixedElectrodePotential3D> fixed_electrodes,
    std::span<const FloatingElectrodeCharge3D> floating_electrodes,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::physics
