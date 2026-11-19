#pragma once

#include "nova/fem/p1_flux_3d.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nova::fem {

// Diagnostics of continuity of the raw, cell-averaged P1 conductive flux.
// For an interior face F, n points out of first_cell and
// jump = (q_first - q_second) dot n. Swapping cell owners leaves jump unchanged.
struct InteriorFluxJumpFace3D {
    BoundaryFace3D nodes{};  // sorted global vertex IDs
    std::size_t first_cell{};
    std::size_t second_cell{};
    double area{};
    double diameter{};      // longest edge of the triangular face
    double normal_jump{};   // q_first dot n_first + q_second dot n_second
};

struct InteriorFluxJumpReport3D {
    std::vector<InteriorFluxJumpFace3D> faces;
    // sqrt(0.5 * sum_{interior F on cell boundary} h_F |F| jump_F^2)
    std::vector<double> cell_indicators;
    // sqrt(sum_{interior F} h_F |F| jump_F^2)
    double global_indicator{};
    double maximum_absolute_normal_jump{};
};

// Each interior face is counted exactly once and its nonnegative squared
// contribution is shared equally between the two incident tetrahedra.
// Boundary faces are excluded. Values measure jumps in averaged flux, NOT
// local conservation, a full residual-based error estimate, or a guaranteed
// upper bound on temperature error. In particular source residuals, boundary
// mismatch and sub-element conductivity variation are not represented.
// The routine validates the mesh and finite cell fluxes. No solve is run.
[[nodiscard]] InteriorFluxJumpReport3D evaluate_interior_flux_jumps_3d(
    const TetrahedralMesh3D& mesh,
    std::span<const P1TetrahedralFlux3D> cell_fields
);

} // namespace nova::fem
