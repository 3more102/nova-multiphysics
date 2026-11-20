#pragma once

#include "nova/fem/transient_heat_3d.hpp"

namespace nova::linalg {
class CachedConjugateGradient;
}

namespace nova::fem {

// Internal overload used by adaptive integration to retain an instance-owned
// IC(0) cache across otherwise independent fixed-step trials.
[[nodiscard]] TransientHeat3DResult solve_transient_heat_3d_with_cache(
    const TransientHeat3DProblem& problem,
    const TransientHeat3DOptions& options,
    nova::linalg::CachedConjugateGradient& cache
);

} // namespace nova::fem
