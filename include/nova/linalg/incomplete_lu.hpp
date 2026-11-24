#pragma once
#include "nova/linalg/csr_matrix.hpp"
namespace nova::linalg {
// Natural-order, zero-fill ILU for general square matrices. No pivoting or
// diagonal shifts. Missing, zero and numerically unresolved pivots are rejected.
class IncompleteLU0 {
public:
    explicit IncompleteLU0(const CsrMatrix& matrix);
    [[nodiscard]] std::vector<double> solve(std::span<const double> rhs) const;
private:
    std::vector<std::size_t> offsets_, columns_, diagonal_;
    std::vector<double> values_;
};
}
