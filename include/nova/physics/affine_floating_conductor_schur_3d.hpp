#pragma once

#include "nova/physics/floating_conductor_schur_3d.hpp"

#include <span>
#include <vector>

namespace nova::physics {

// P1 floating-conductor Schur solver including impressed volume charge and
// prescribed normal electric displacement on selected non-Dirichlet faces.
//
// The homogeneous terminal matrix C is unchanged by the linear sources.
// The affine electrode reactions are Q = C*V + Q0, where Q0 is calculated
// from a loaded, zero-terminal-voltage FEM solve as K*V0 - F at electrode
// Dirichlet nodes (using the original diffusion solver's quadrature).
//
// All non-electrode potential constraints become grounded references, as
// in solve_floating_conductor_schur_3d. No Robin/nonlinear dielectric,
// conductor overlap, or solution certification is implied.
struct AffineSchurFloatingConductor3DResult {
    CapacitanceMatrix3D discrete_capacitance;
    FloatingElectrodePotentials3D terminals;
    // Sorted with terminals.electrode_ids. Q0 in the affine relation.
    std::vector<double> impressed_charge_offset_coulombs;
    Electrostatics3DResult field;
    // Loaded weak-form reaction K*V - F, NOT merely K*V.
    std::vector<DiscreteElectrodeReaction3D> loaded_discrete_reactions;
    // Raw boundary-D diagnostic ignores FEM volume/source load accounting.
    std::vector<ElectrodeCharge3D> boundary_flux_charges;
    double max_floating_reaction_mismatch_coulombs{};
    double max_flux_reaction_disagreement_coulombs{};
};

[[nodiscard]] AffineSchurFloatingConductor3DResult
solve_affine_floating_conductor_schur_3d(
    const Electrostatics3DProblem& dielectric_problem,
    std::span<const ElectrodeBoundaryGroup3D> electrode_groups,
    std::span<const FixedElectrodePotential3D> fixed_electrodes,
    std::span<const FloatingElectrodeCharge3D> floating_electrodes,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::physics
