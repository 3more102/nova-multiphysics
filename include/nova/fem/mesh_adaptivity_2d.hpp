#pragma once

#include "nova/fem/h1_space_2d.hpp"
#include "nova/fem/error_estimator_2d.hpp"
#include "nova/fem/mesh_refinement_2d.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace nova::fem {

// Elementwise Zienkiewicz-Zhu-style gradient-recovery indicator. The
// recovered nodal gradient is an area-weighted average of element-centroid
// gradients. This is a refinement heuristic, not a guaranteed error bound.
[[nodiscard]] std::vector<double> estimate_h1_gradient_recovery_indicators_2d(
    const TriangleMesh2D& mesh,
    H1TriangleOrder order,
    std::span<const double> dof_values
);

// Dörfler bulk marking: returns the smallest deterministic set, ordered by
// descending indicator, whose squared indicators reach theta times the
// global squared-indicator sum. theta must lie in (0,1].
[[nodiscard]] std::vector<std::size_t> mark_elements_dorfler_2d(
    std::span<const double> indicators,
    double theta
);

// Selective conforming refinement. Each marked element requests its longest
// edge; all cells incident to a requested edge are split, preventing hanging
// nodes. Cells with one, two, or three split edges create two, three, or four
// children respectively. Boundary edges are split only when their edge is
// selected, preserving the original directed traversal.
[[nodiscard]] RefinedTriangleMesh2D refine_triangles_by_marked_elements_2d(
    const TriangleMesh2D& mesh,
    const std::vector<BoundaryEdge2D>& boundary_edges,
    std::span<const std::size_t> marked_elements
);


struct P1ResidualAdaptation2D {
    DiffusionResidualIndicators2D estimate;
    std::vector<std::size_t> marked_elements;
    RefinedTriangleMesh2D refinement;
};

// Runs one residual-driven adaptation cycle: estimate, Dörfler-mark, then
// conformingly refine the marked longest edges. An all-zero estimator returns
// an unchanged mesh with an empty marked-element list.
[[nodiscard]] P1ResidualAdaptation2D adapt_p1_diffusion_mesh_by_residual_2d(
    const TriangleMesh2D& mesh,
    const std::vector<BoundaryEdge2D>& boundary_edges,
    std::span<const double> nodal_values,
    std::span<const double> element_diffusivity,
    const std::function<double(double, double)>& source,
    double theta
);

} // namespace nova::fem
