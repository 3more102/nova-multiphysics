#pragma once

#include "nova/fem/transient_heat_2d.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace nova::fem::detail {

struct TransientHeat2DBoundaryData {
    std::vector<Point2D> coordinates;
    std::vector<std::vector<TransientScalarField2D>> functions;
    std::vector<bool> constrained;
};

struct PreparedTransientHeat2D {
    ScalarH1Operators2D operators;
    TransientHeat2DBoundaryData boundary;
    std::vector<double> capacity_weights;
};

// Each adaptive trial slot owns its dt-dependent base operator and CG cache.
// The Robin contribution is time-dependent and is added at each time level;
// CachedConjugateGradient reuses IC(0) only for an exactly identical system.
struct TransientHeat2DStepCache {
    std::optional<linalg::CsrMatrix> base_system;
    double time_step{};
    TransientTimeScheme scheme{TransientTimeScheme::BackwardEuler};
    linalg::CachedConjugateGradient cached_cg;
    std::size_t base_operator_builds{};
};

[[nodiscard]] PreparedTransientHeat2D prepare_transient_heat_2d(
    const TransientHeat2DProblem& problem,
    const H1FunctionSpace2D& space
);

[[nodiscard]] std::vector<double> initial_heat_state(
    const TransientHeat2DProblem& problem,
    const H1FunctionSpace2D& space
);

[[nodiscard]] TransientHeat2DResult solve_prepared_transient_heat_2d(
    const TransientHeat2DProblem& problem,
    const H1FunctionSpace2D& space,
    const PreparedTransientHeat2D& prepared,
    const TransientHeat2DOptions& options,
    std::vector<double> initial_state,
    TransientHeat2DStepCache& cache
);

} // namespace nova::fem::detail
