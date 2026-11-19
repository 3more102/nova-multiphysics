#include "nova/fem/diffusion_2d.hpp"

namespace nova::fem {

// Keep the established P1 API while routing all scalar 2D diffusion
// assembly through the shared H1 P1/P2 kernel.
ScalarDiffusion2DResult solve_scalar_diffusion_2d(
    const ScalarDiffusion2DProblem& problem,
    const linalg::SolverOptions& solver_options
) {
    return solve_scalar_diffusion_h1_2d(
        problem, H1TriangleOrder::P1, solver_options
    );
}

} // namespace nova::fem
