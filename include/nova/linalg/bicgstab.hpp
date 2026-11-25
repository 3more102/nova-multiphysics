#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/solver_types.hpp"

#include <span>

namespace nova::linalg {

// Solves A*x=b using BiCGSTAB for general nonsymmetric square systems.
// A successful solve satisfies the recomputed physical residual criterion.
// Krylov corrections use power-of-two residual scaling and angle-based
// breakdown tests. Matvecs, coefficients and the solution must fit in double.
// Jacobi requires a finite, nonzero diagonal with representable reciprocal.
[[nodiscard]] SolverResult solve_bicgstab(
    const CsrMatrix& matrix,
    std::span<const double> right_hand_side,
    const SolverOptions& options = {},
    std::span<const double> initial_guess = {}
);

} // namespace nova::linalg
