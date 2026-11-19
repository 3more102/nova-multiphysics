#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include "nova/linalg/incomplete_cholesky.hpp"
#include "nova/linalg/solver_types.hpp"

#include <optional>
#include <span>

namespace nova::linalg {

// Solves A*x=b using conjugate gradients. A must be symmetric positive-definite.
// Jacobi preconditioning requires a finite, strictly positive matrix diagonal.
// IncompleteCholesky0 builds one unshifted, zero-fill factor per solve.
[[nodiscard]] SolverResult solve_conjugate_gradient(
    const CsrMatrix& matrix,
    std::span<const double> right_hand_side,
    const SolverOptions& options = {},
    std::span<const double> initial_guess = {}
);


// Keeps an IC(0) factor for repeated solves against an exactly identical CSR
// matrix. Any changed dimension, sparsity pattern or coefficient rebuilds it.
// Options and right-hand sides are still validated separately on each solve.
// Own one cache per simulation; no global mutable solver state is used.
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
