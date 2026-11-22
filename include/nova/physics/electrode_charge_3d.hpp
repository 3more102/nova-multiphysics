#pragma once

#include "nova/physics/electrostatics_3d.hpp"
#include "nova/linalg/solver_types.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nova::physics {

// Approximate charge on a prescribed-potential conductor, using the P1
// element-averaged electric displacement on the adjacent dielectric side.
// Positive conductor charge points into the dielectric: Q = -integral(D.n) dA,
// where n is the outward normal of the DIELECTRIC computational domain.
// This raw boundary-flux diagnostic is not locally equilibrated or guaranteed
// to obey Gauss's law on coarse/unstructured meshes.
struct ElectrodeCharge3D {
    int electrode_id{};
    std::size_t face_count{};
    double area_square_metres{};
    double potential_volts{};
    double dielectric_outward_charge_coulombs{};
    double conductor_charge_coulombs{};
};

// Faces must be unique exterior triangles with vertices all fixed by the
// problem's prescribed-potential conditions. Every face belonging to one
// electrode ID must have a common potential to floating-point precision.
// The returned groups are sorted by nonnegative electrode ID. Winding is
// ignored; geometric normals are oriented outward from the dielectric.
// Neither the fields nor the mesh are re-solved or modified.
[[nodiscard]] std::vector<ElectrodeCharge3D> evaluate_electrode_charges_3d(
    const Electrostatics3DProblem& problem,
    const Electrostatics3DResult& solution,
    std::span<const fem::BoundaryFace3D> electrode_faces,
    std::span<const int> electrode_ids
);


struct ElectrodeBoundaryGroup3D {
    int electrode_id{};
    std::vector<fem::BoundaryFace3D> faces;
};

struct CapacitanceMatrix3D {
    // Sorted electrode IDs. matrix_farads is row-major: C(i,j) is the
    // conductor charge on electrode i for a unit voltage on electrode j.
    std::vector<int> electrode_ids;
    std::vector<double> matrix_farads;
    // Raw-flux reciprocity diagnostics; these are not symmetrized.
    double max_reciprocity_defect_farads{};
    double max_relative_reciprocity_defect{};
};

// Extract the experimental Maxwell capacitance matrix using one unit-voltage
// solve per electrode. Existing non-electrode potential boundaries become
// grounded references. Free-charge and prescribed-displacement loads are
// excluded from this linear capacitance calculation.
[[nodiscard]] CapacitanceMatrix3D extract_capacitance_matrix_3d(
    const Electrostatics3DProblem& dielectric_problem,
    std::span<const ElectrodeBoundaryGroup3D> electrodes,
    const linalg::SolverOptions& solver_options = {}
);


// Discrete weak-form reaction charge for equipotential Dirichlet electrodes.
// For source-free 3D P1 electrostatics, the net conductor reaction is
//   Q_g = sum_{i in electrode g} (K V)_i
//       = sum_K integral_K epsilon grad(w_g) dot grad(V) dV,
// where w_g is one on the electrode's selected boundary vertices and zero
// elsewhere. This is a Galerkin reaction, distinct from direct face-D flux.
// Material quadrature matches the underlying P1 FEM assembly.
//
// The selected triangle groups must cover each electrode's entire contact;
// no selected node may belong to multiple electrode IDs. Every selected
// vertex must be explicitly fixed by the source problem. The method rejects
// free-charge callbacks, cell volume charge, and prescribed displacement
// loads: including them requires reaction load accounting.
//
// This is an algebraic FEM consistency diagnostic, not an independently
// converged physical surface-charge measurement.
struct DiscreteElectrodeReaction3D {
    int electrode_id{};
    std::size_t selected_face_count{};
    double potential_volts{};
    double reaction_charge_coulombs{};
};

[[nodiscard]] std::vector<DiscreteElectrodeReaction3D>
evaluate_discrete_electrode_reactions_3d(
    const Electrostatics3DProblem& problem,
    const Electrostatics3DResult& solution,
    std::span<const fem::BoundaryFace3D> electrode_faces,
    std::span<const int> electrode_ids
);


// Re-extract a capacitance matrix using the discrete P1 weak-form
// conductor reaction (K*V)_electrode rather than the element-averaged D.n
// surface-flux integral used by extract_capacitance_matrix_3d.
// The same validated electrode sets, charge-free unit-voltage basis solves,
// and optional grounded non-electrode references are used; it performs
// ONE solve per electrode, not an additional solve for each diagnostic.
// Symmetry follows from the Galerkin bilinear form up to linear-solver and
// quadrature roundoff. No hidden symmetrization or gauge fixing is applied.
// This is a discretization-consistent capacitance approximation, not a
// certified continuum Maxwell capacitance or a floating-boundary FEM solver.
[[nodiscard]] CapacitanceMatrix3D extract_discrete_capacitance_matrix_3d(
    const Electrostatics3DProblem& dielectric_problem,
    std::span<const ElectrodeBoundaryGroup3D> electrodes,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::physics
