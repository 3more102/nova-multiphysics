#pragma once

#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace nova::linalg {

// Natural-order, zero-fill incomplete Cholesky for a symmetric matrix.
// The lower factor has exactly the input's lower-triangular sparsity pattern.
// Missing diagonal entries, asymmetric patterns/values, and unresolved
// pivots are rejected. No reordering, diagonal shift, or fallback is applied.
class IncompleteCholesky0 {
public:
    explicit IncompleteCholesky0(const CsrMatrix& matrix);

    [[nodiscard]] std::size_t dimension() const noexcept { return diagonal_.size(); }
    [[nodiscard]] std::size_t stored_entries() const noexcept {
        return lower_.size() + dimension();
    }

    // Apply (L L^T)^-1 to a finite right-hand side. A factor is reusable.
    [[nodiscard]] std::vector<double> apply(std::span<const double> rhs) const;

private:
    std::vector<std::size_t> row_offsets_;
    std::vector<std::size_t> columns_;
    std::vector<double> lower_;
    std::vector<double> diagonal_;
};

} // namespace nova::linalg
