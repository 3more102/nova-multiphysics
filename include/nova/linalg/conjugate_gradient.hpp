#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/solver_types.hpp"

#include <span>

namespace nova::linalg {

// Solves A*x=b using conjugate gradients. A must be symmetric positive-definite.
// Jacobi preconditioning requires a finite, strictly positive matrix diagonal.
// IC(0) builds an unshifted natural-order zero-fill factor for each CG solve.
[[nodiscard]] SolverResult solve_conjugate_gradient(
    const CsrMatrix& matrix,
    std::span<const double> right_hand_side,
    const SolverOptions& options = {},
    std::span<const double> initial_guess = {}
);

} // namespace nova::linalg
