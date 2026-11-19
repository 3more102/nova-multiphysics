#pragma once

#include "nova/fem/diffusion_3d.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace nova::fem {

// Signed face-integrated heat rate [W], positive outward from owner_cell.
// A shared interior face appears exactly once: the neighboring cell uses
// the negative of outward_heat_rate. Vertices are in sorted key order.
enum class FaceFlowBoundaryKind3D {
    interior, dirichlet, neumann, robin, insulated
};

struct ConservativeFaceFlow3D {
    BoundaryFace3D nodes{};
    std::size_t owner_cell{};
    std::optional<std::size_t> neighbor_cell{};
    FaceFlowBoundaryKind3D kind{FaceFlowBoundaryKind3D::interior};
    double raw_outward_heat_rate{};
    double outward_heat_rate{};
};

struct ConservativeHeatFlow3D {
    std::vector<ConservativeFaceFlow3D> faces;
    std::vector<double> generated_heat_per_cell;
    // Generated heat minus the sum of outward face rates for each tetrahedron.
    std::vector<double> defect_before;
    std::vector<double> defect_after;
    double maximum_defect_before{};
    double maximum_defect_after{};
    double maximum_absolute_face_correction{};
    double total_generated_heat{};
    double total_exterior_outward_heat{};
};

// Algebraically conservative face-flow projection for a P1 tetrahedral field.
//
// Starts from the two-cell average of raw -K grad(T) on interior faces.
// Preserves the triangle-quadrature-integrated Neumann and Robin rates and
// leaves unassigned exterior faces insulated. Corrects interior-face rates
// along a deterministic spanning forest, placing any remaining per-component
// heat imbalance on one explicitly marked Dirichlet exterior face.
// Components without such an outlet must have compatible source/boundary loads.
//
// Conservation, single-valued interface rates and prescribed boundary rates
// hold to floating-point accuracy. This is NOT an RT0/H(div) minimizer,
// pointwise physical flux, unique distribution of Dirichlet reactions, or a
// substitute for the independent FEM weak-residual convergence audit.
// The correction is tree-dependent; inspect maximum_absolute_face_correction.
//
// Validates the same problem/field contracts as evaluate_diffusion_energy_balance_3d.
// Does not modify or solve the original FEM system.
[[nodiscard]] ConservativeHeatFlow3D reconstruct_conservative_heat_flow_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> nodal_values
);

} // namespace nova::fem
