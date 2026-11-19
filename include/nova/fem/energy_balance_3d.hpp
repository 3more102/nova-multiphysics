#pragma once

#include "nova/fem/diffusion_3d.hpp"

#include <cstddef>
#include <span>

namespace nova::fem {

// Global balance of the *assembled weak form*, evaluated with the same
// P1 tetrahedral/triangular quadrature as solve_scalar_diffusion_3d.
// SI units are W when coordinates and input fields use SI.
// Positive outward rates mean heat leaves the modeled solids.
struct DiffusionEnergyBalance3D {
    double generated_heat{};                 // Integral of volumetric source, signed
    double neumann_outward_heat{};           // Prescribed exterior Neumann heat flow
    double robin_outward_heat{};             // Integral h * (T_h - ambient)
    double constraint_outward_heat{};        // Negative sum of constrained-node FE residuals
    double signed_balance_defect{};          // generated - Neumann - Robin - constraints
    double max_free_equation_residual{};     // Maximum absolute unconstrained weak residual [W]
    std::size_t constrained_nodes{};
    std::size_t free_nodes{};
};

// Reassembles the local FEM weak residual, not the raw element-face flux.
// The constrained-node residual is the algebraic reaction of ALL prescribed
// vertex degrees of freedom, including explicitly fixed interior vertices.
// It is not an individually resolved physical flux on each Dirichlet face.
// An arbitrary trial field may have a nonzero signed defect and large free
// residuals; no convergence claim is made by this function.
// Validates any supplied material/source callbacks, face assignments,
// prescribed data, finite nodal values, and fixed-node field consistency.
// Does not modify the input, invoke a solver, or claim local conservation.
[[nodiscard]] DiffusionEnergyBalance3D evaluate_diffusion_energy_balance_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> nodal_values
);

} // namespace nova::fem
