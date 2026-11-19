#pragma once

#include "nova/linalg/csr_matrix.hpp"

namespace nova::linalg {

// Natural-order IC(0) for a symmetric matrix. The lower factor has exactly
// the input's lower-triangular storage pattern, including its diagonal.
// A missing diagonal, asymmetric pattern/value or unresolved pivot raises
// std::invalid_argument. No ordering, diagonal shift or fallback is applied.
class IncompleteCholesky0 {
public:
    explicit IncompleteCholesky0(const CsrMatrix& matrix);

    [[nodiscard]] std::size_t dimension() const noexcept { return diagonal_.size(); }
    [[nodiscard]] std::size_t stored_entries() const noexcept { return lower_.size() + dimension(); }

    // Solve L L^T z = r. The factor can be reused for multiple right-hand sides.
    [[nodiscard]] std::vector<double> apply(std::span<const double> rhs) const;

private:
    std::vector<std::size_t> row_offsets_;
    std::vector<std::size_t> columns_;
    std::vector<double> lower_;
    std::vector<double> diagonal_;
};

} // namespace nova::linalg
