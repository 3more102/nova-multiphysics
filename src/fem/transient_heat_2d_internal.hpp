#pragma once

#include "nova/fem/transient_heat_2d.hpp"

#include <optional>
#include <vector>

namespace nova::fem::detail {

struct BoundaryData {
    std::vector<Point2D> coordinates;
    std::vector<std::vector<TransientScalarField2D>> functions;
    std::vector<bool> constrained;
};

struct PreparedTransientHeat2D {
    ScalarH1Operators2D operators;
    BoundaryData boundary;
    std::vector<double> capacity_weights;
};

// One cache belongs to either the full-step or half-step adaptive trial slot.
// The base M + theta*dt*K operator is retained for an identical time step.
// The exact constrained matrix, including the current Robin matrix, remains
// the authority for IC(0) reuse inside CachedConjugateGradient.
struct TransientHeat2DStepCache {
    std::optional<linalg::CsrMatrix> base_system;
    double time_step{};
    TransientTimeScheme scheme{TransientTimeScheme::BackwardEuler};
    linalg::CachedConjugateGradient cached_cg;
    std::size_t operator_builds{};
};

[[nodiscard]] PreparedTransientHeat2D prepare_transient_heat_2d(
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
