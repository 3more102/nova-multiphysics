#pragma once

#include "nova/fem/mesh_refinement_2d.hpp"
#include "nova/fem/transient_heat_2d.hpp"

#include <span>
#include <vector>

namespace nova::fem {

// Transfers a continuous P1/P2 H1 field to a conforming midpoint refinement
// (uniform or red-green adaptive). P2 transfer evaluates the original
// quadratic element polynomial at every new fine-edge midpoint.
[[nodiscard]] std::vector<double> prolongate_h1_refined_2d(
    const TriangleMesh2D& coarse_mesh,
    const RefinedTriangleMesh2D& refined,
    H1TriangleOrder order,
    std::span<const double> coarse_dof_values
);

[[nodiscard]] std::vector<double> prolongate_h1_uniform_2d(
    const TriangleMesh2D& coarse_mesh,
    const RefinedTriangleMesh2D& refined,
    H1TriangleOrder order,
    std::span<const double> coarse_dof_values
);

struct RefinedTransientHeat2DProblem {
    TransientHeat2DProblem problem;
    RefinedTriangleMesh2D refinement;
};

// Uniformly refines a transient heat problem, splits every assigned boundary
// edge while preserving its direction and callback, and transfers an optional
// exact P1/P2 restart vector. Spatial and time-dependent field callbacks are
// retained unchanged because they are evaluated in physical coordinates.
[[nodiscard]] RefinedTransientHeat2DProblem
refine_transient_heat_problem_uniformly_2d(
    const TransientHeat2DProblem& problem,
    H1TriangleOrder order
);

// Applies selective red-green refinement to a transient heat problem. Boundary
// edges are split only when their parent edge is bisected; other edges and
// callbacks are preserved. Optional restart states are transferred for P1/P2.
[[nodiscard]] RefinedTransientHeat2DProblem
refine_transient_heat_problem_adaptively_2d(
    const TransientHeat2DProblem& problem,
    H1TriangleOrder order,
    std::span<const std::size_t> marked_elements
);

} // namespace nova::fem
