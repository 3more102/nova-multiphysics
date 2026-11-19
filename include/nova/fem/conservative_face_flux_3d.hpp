#pragma once

#include "nova/fem/p1_flux_3d.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace nova::fem {

// Face-integrated rate: positive outward from first_cell. The optional
// neighbor receives precisely the negative of this shared interior rate.
struct ConservativeFaceHeatRate3D {
    BoundaryFace3D vertices{}; // sorted global node IDs
    std::size_t first_cell{};
    std::optional<std::size_t> second_cell{};
    double area{};
    double outward_heat_rate{};
};

struct PrescribedBoundaryHeatRate3D {
    BoundaryFace3D face{};
    double outward_heat_rate{}; // positive outward [W]
};

struct ConservativeFaceFlux3DResult {
    std::vector<ConservativeFaceHeatRate3D> faces;
    std::vector<double> cell_balance_defect; // source - sum(outward) [W]
    double maximum_absolute_cell_defect{};
    std::size_t correction_iterations{};
};

// Sparse graph-Laplacian projection of face-integrated rates from discontinuous
// raw P1 flux. Interior faces have a single equal/opposite heat-flow unknown;
// unassigned exterior faces are insulated. Adjustable boundary triangles are
// the only unspecified exterior outlets (e.g., physical Dirichlet boundaries).
// Incompatible all-Neumann component budgets are rejected.
// Rates are conservative to the documented numerical tolerance, not a
// reconstructed pointwise H(div) field or an a posteriori error certificate.
[[nodiscard]] ConservativeFaceFlux3DResult reconstruct_conservative_face_flux_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const P1TetrahedralFlux3D> raw_cell_flux,
    std::span<const double> cell_source_heat_rate,
    std::span<const PrescribedBoundaryHeatRate3D> prescribed_boundary_rates,
    std::span<const BoundaryFace3D> adjustable_boundary_faces
);

// Integrates source, Neumann and Robin data using solver-matching degree-two
// quadrature. Fully Dirichlet exterior triangles become adjustable outlets;
// all other unassigned boundary faces are insulated.
[[nodiscard]] ConservativeFaceFlux3DResult reconstruct_conservative_face_flux_3d(
    const ScalarDiffusion3DProblem& problem,
    std::span<const double> nodal_values
);

} // namespace nova::fem
