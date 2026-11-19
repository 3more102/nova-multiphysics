#pragma once

#include "nova/fem/diffusion_3d.hpp"

#include <span>

namespace nova::fem {

// Globally integrated steady energy budget in watts when using SI units.
// The Dirichlet term is the NEGATIVE assembled weak-form reaction at physical
// prescribed nodes. It is not a pointwise reconstructed boundary heat flux.
// All terms use the same quadrature rules and material choices as the FEM
// assembly. A converged free-node solution closes the discrete balance up to
// numerical error; this is not a local H(div)-equilibrated flux certificate.
struct SteadyEnergyBalance3D {
    double volumetric_generation{};
    double neumann_outward{};
    double robin_outward{};
    double dirichlet_outward{};
    double total_outward{};
    // total_outward - volumetric_generation; ideally zero.
    double imbalance{};
};

// Independently re-evaluates the assembled weak form for a supplied P1 field.
// Explicit Dirichlet nodes and vertices of Dirichlet faces count as physical
// constraints; artificial pure-Neumann solver gauge pins do NOT.
// All boundary faces must be distinct valid exterior faces. Nonfinite inputs,
// missing material/source data and inconsistent vector dimensions are rejected.
[[nodiscard]] SteadyEnergyBalance3D evaluate_steady_energy_balance_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> temperatures
);

} // namespace nova::fem
