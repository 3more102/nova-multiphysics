#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/mesh_refinement_3d.hpp"
#include "nova/fem/residual_error_estimator_3d.hpp"

#include <vector>

namespace nova::fem {

// One explicit solve -> estimate -> Dörfler mark -> conforming local refine
// -> transfer -> resolve cycle for stationary P1 tetrahedral diffusion.
//
// Restricted to per-cell positive scalar conductivity and a spatial source
// callback (or zero source). Dirichlet vertex constraints and face callbacks,
// prescribed Neumann faces and otherwise insulated exterior faces are
// transferred. Cell sources, tensor/spatial conductivity and Robin data
// are intentionally rejected until their estimators/transfer are supported.
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

// Bounded sequence of complete solve/estimate/mark/refine/resolve steps.
// Each subsequent step consumes the previous step's refined PDE/mesh.
// The full history is retained for inspection. No estimator tolerance,
// mesh-quality guarantee, or certified error bound is implied.
struct AdaptiveDiffusionCycles3D {
    std::vector<AdaptiveDiffusionStep3D> steps;
    bool stopped_without_marking{};
};

// Executes at most max_cycles steps (max_cycles must be >= 1), stopping
// early if a step has no marked tetrahedra. This deliberately recomputes
// each step's coarse solve; matrix/iterate reuse is a future optimization.
// The final field and mesh are steps.back().refined_solution and
// steps.back().refined_problem.
[[nodiscard]] AdaptiveDiffusionCycles3D run_adaptive_diffusion_cycles_3d(
    const ScalarDiffusion3DProblem& problem,
    double dorfler_theta,
    std::size_t max_cycles,
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::fem
