#pragma once

#include "nova/linalg/csr_matrix.hpp"
#include <cstddef>
#include <vector>

namespace nova::linalg {

struct SymmetricEigenOptions {
    std::size_t max_dimension{512}; // deliberately bounded dense reference solver
    std::size_t max_sweeps{80};
    double relative_tolerance{1e-13};
};

struct SymmetricEigenpair {
    double eigenvalue{};
    std::vector<double> eigenvector; // M-normalized
    double relative_residual{}; // ||Kx-lambda Mx|| / (||K||_F+|lambda| ||M||_F)||x||
};

struct SymmetricEigenResult {
    std::vector<SymmetricEigenpair> eigenpairs; // ascending eigenvalues
    std::size_t sweeps{};
    bool converged{};
};

// Real symmetric K, symmetric positive definite M. Cholesky reduction
// followed by cyclic Jacobi rotations. O(n^3) time, O(n^2) storage.
// Throws for invalid/nonfinite inputs, non-SPD mass or dimension limit.
[[nodiscard]] SymmetricEigenResult solve_symmetric_generalized_eigen(
    const CsrMatrix& stiffness,const CsrMatrix& mass,
    const SymmetricEigenOptions& options={});

} // namespace nova::linalg
