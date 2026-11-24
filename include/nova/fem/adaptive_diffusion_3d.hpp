#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/mesh_refinement_3d.hpp"
#include "nova/fem/residual_error_estimator_3d.hpp"

#include <vector>

namespace nova::fem {

// One explicit solve -> estimate -> Dörfler mark -> conforming local refine
// -> transfer -> resolve cycle for stationary P1 tetrahedral diffusion.
//
// Supports any one of the scalar/tensor spatial or cellwise conductivity
// representations accepted by ScalarDiffusion3DProblem, with a spatial source
// callback (or zero source), with optional additive cellwise source vectors.
// Dirichlet, Neumann, Robin and insulated exterior faces are transferred.
// Flux-jump indicators use the full reconstructed physical P1 cell flux.
//
// A zero indicator field causes no refinement and no redundant second solve.
// This is an experimental error-indicator workflow, NOT an error bound.
struct AdaptiveDiffusionStep3D {
    ScalarDiffusion3DResult coarse_solution;
    DiffusionResidualIndicators3D coarse_indicators;
    IndicatorRefinementResult3D adaptation;
    std::vector<double> transferred_p1_values;
    ScalarDiffusion3DProblem refined_problem;
    ScalarDiffusion3DResult refined_solution;
    DiffusionResidualIndicators3D refined_indicators;
};

[[nodiscard]] AdaptiveDiffusionStep3D solve_adaptive_diffusion_step_3d(
    const ScalarDiffusion3DProblem& problem,
    double dorfler_theta,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::fem
