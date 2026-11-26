#pragma once

#include "nova/physics/floating_electrodes_3d.hpp"

#include <span>
#include <vector>

namespace nova::physics {

// Experimental P1 Schur-complement elimination for LINEAR floating
// conductors, including prescribed volume free charge and Neumann
// displacement loads. Unit-electrode Dirichlet FEM solves implicitly
// eliminate free interior nodes to form C = K_ee - K_ei K_ii^{-1} K_ie.
// A loaded all-zero-terminal solve gives the affine charge offset Q0.
// Prescribed floating charges solve C_FF V_F = Q_F - Q0_F - C_FK V_K.
// A final FEM solve at those potentials returns spatial E and D fields.
//
// This is a discrete Schur complement assembled through P1 basis solves,
// not an explicitly assembled mixed/saddle-point solver. Ground references
// external to electrode groups are set to zero volts, matching the
// capacitance extraction API. Impressed volume/Neumann loads contribute to
// Q0 via a matching discrete K*V-rhs reaction diagnostic. Neumann faces
// cannot overlap Dirichlet electrode faces. For loaded cases, the usual
// source-free sum(V*Q)=2*U energy identity no longer applies.
// A floating conductor with no voltage reference is rejected if its
// Schur block has a gauge nullspace.
struct SchurFloatingConductor3DResult {
    CapacitanceMatrix3D discrete_capacitance;
    // Loaded zero-voltage terminal reactions Q0, sorted by electrode ID [C].
    // All entries are zero when no impressed charge or Neumann loads exist.
    std::vector<double> load_offset_charge_coulombs;
    // Estimated terminal charges include the affine Q0 offset.
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
