#pragma once

#include "nova/fem/diffusion_3d.hpp"

#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace nova::fem {

// A unique triangular face, listed once in canonical vertex order.
// Heat rate is positive outward from owner_cell; on an interior face the
// neighbor's outward heat rate is its negative (exactly one shared unknown).
// Boundary faces have neighbor_cell == no_neighbor.
struct ConservativeFaceHeat3D {
    static constexpr std::size_t no_neighbor = std::numeric_limits<std::size_t>::max();
    BoundaryFace3D vertices{};
    std::size_t owner_cell{};
    std::size_t neighbor_cell{no_neighbor};
    double outward_heat_rate{}; // integrated W, NOT a pointwise flux density
};

struct ConservativeHeatFluxResult3D {
    std::vector<ConservativeFaceHeat3D> faces;
    std::vector<double> generated_heat_per_cell;
    std::vector<double> outward_heat_per_cell;
    double max_cell_balance_defect{};  // max |generation - sum(outward faces)| [W]
    std::size_t correction_iterations{};
};

// Reconstruct continuous, cell-equilibrated face-integrated heat rates from
// a STEADY P1 temperature field. Start from averaged cell fluxes; minimize the
// unweighted sum of squared changes to free face rates subject to one heat
// balance equation per tetrahedron. On interior faces, one rate is shared;
// Neumann and Robin rates are integrated and held fixed; unassigned exterior
// faces are insulated (fixed zero). Explicit Dirichlet FACES are free outlets.
//
// Uses the same source and boundary quadrature as the stationary FEM solver.
// Pure-Neumann connected components must have balanced prescribed heat;
// otherwise the constrained face-flow problem is infeasible and throws.
// Vertex-only Dirichlet constraints, including interior pins, are NOT inferred
// as physical boundary outlets. Explicitly label exterior Dirichlet faces
// when heat needs to leave through them.
//
// This is a face-INTEGRATED, topological conservative postprocessor, NOT an
// H(div) vector-field or RT0 reconstruction and NOT a transient energy audit.
// No changes are made to the input problem or temperature field.
[[nodiscard]] ConservativeHeatFluxResult3D reconstruct_conservative_heat_flux_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> nodal_values
);

} // namespace nova::fem
