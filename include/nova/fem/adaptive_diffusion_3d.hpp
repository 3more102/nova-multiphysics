#pragma once

#include "nova/fem/diffusion_3d.hpp"
#include "nova/fem/mesh_refinement_3d.hpp"
#include "nova/fem/residual_error_estimator_3d.hpp"

#include <cstddef>
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

// Controls an automatic sequence of one-cycle adaptive solves. Tolerances are
// applied to the residual indicator, not to a certified continuum error.
struct AdaptiveDiffusionOptions3D {
    double dorfler_theta{0.5};
    double absolute_tolerance{1.0e-8};
    double relative_tolerance{1.0e-3};
    std::size_t max_cycles{8};
};

enum class AdaptiveDiffusionStopReason3D {
    ToleranceSatisfied,
    NoRefinementRequired,
    MaximumCyclesReached
};

struct AdaptiveDiffusionCycleSummary3D {
    std::size_t cycle{};
    std::size_t coarse_nodes{};
    std::size_t coarse_elements{};
    std::size_t marked_elements{};
    std::size_t refined_nodes{};
    std::size_t refined_elements{};
    double coarse_indicator{};
    double refined_indicator{};
};

struct AdaptiveDiffusionRun3D {
    ScalarDiffusion3DProblem final_problem;
    ScalarDiffusion3DResult final_solution;
    DiffusionResidualIndicators3D final_indicators;
    std::vector<AdaptiveDiffusionCycleSummary3D> cycles;
    double initial_indicator{};
    double target_indicator{};
    bool tolerance_satisfied{};
    AdaptiveDiffusionStopReason3D stop_reason{
        AdaptiveDiffusionStopReason3D::MaximumCyclesReached};
};

// Runs up to max_cycles, stopping when the final indicator is no larger than
// max(absolute_tolerance, relative_tolerance * initial_indicator), or when no
// further elements can be marked. The initial indicator is the first coarse
// solve's total_error; the final result contains the last refined state.
[[nodiscard]] AdaptiveDiffusionRun3D solve_adaptive_diffusion_3d(
    const ScalarDiffusion3DProblem& problem,
    const AdaptiveDiffusionOptions3D& options = {},
    const linalg::SolverOptions& solver_options = {}
);

} // namespace nova::fem
