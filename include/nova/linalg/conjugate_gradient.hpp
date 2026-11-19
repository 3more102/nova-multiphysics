#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/incomplete_cholesky.hpp"
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

// Instance-owned cache for repeated CG solves. IC(0) is reused only when
// dimensions, sparsity and every stored coefficient are exactly identical.
// Any changed matrix rebuilds the factor; non-IC(0) options use ordinary CG.
// Keep one cache per simulation/integration; instances are not thread-safe.
class CachedConjugateGradient {
public:
    [[nodiscard]] SolverResult solve(
        const CsrMatrix& matrix,
        std::span<const double> right_hand_side,
        const SolverOptions& options = {},
        std::span<const double> initial_guess = {}
    );

    [[nodiscard]] std::size_t factor_builds() const noexcept { return factor_builds_; }

private:
    std::optional<CsrMatrix> cached_matrix_;
    std::optional<IncompleteCholesky0> cached_factor_;
    std::size_t factor_builds_{0};
};

} // namespace nova::linalg
