#pragma once

#include "nova/linalg/csr_matrix.hpp"

#include <cstddef>
#include <vector>

namespace nova::linalg {

struct SparseSymmetricEigenOptions {
    std::size_t max_modes{6};
    std::size_t max_subspace{64}; // dense Rayleigh-Ritz projection is bounded to 512
    std::size_t max_mass_iterations{2000};
    double relative_tolerance{1.0e-8};
    double mass_solve_relative_tolerance{1.0e-12};
    double breakdown_tolerance{1.0e-12};
};

struct SparseSymmetricEigenpair {
    double eigenvalue{};
    std::vector<double> eigenvector; // normalized so x^T M x = 1
    double relative_residual{}; // ||Kx-lambda Mx|| / ((||K||_F+|lambda|||M||_F)||x||)
};

struct SparseSymmetricEigenResult {
    std::vector<SparseSymmetricEigenpair> eigenpairs; // descending algebraic eigenvalue
    std::size_t subspace_dimension{};
    std::size_t operator_applications{};
    std::size_t mass_solve_iterations{};
    bool converged{};
};

// Computes the largest algebraic eigenvalues of K*x=lambda*M*x using a
// fully reorthogonalized, generalized Krylov/Rayleigh-Ritz method.
//
// K must be real symmetric and M symmetric positive definite. M systems are
// solved by the existing Jacobi-preconditioned conjugate-gradient solver.
// The projected dense problem is bounded to 512 dimensions; the full sparse
// matrices and eigenvectors are not converted to dense storage.
//
// A deterministic block of initial vectors helps resolve repeated eigenvalues.
// Exhausted invariant subspaces are supplemented by orthogonal coordinate seeds
// while projection capacity remains. Residual convergence alone cannot certify
// global spectral ordering for a projection smaller than the full dimension.
// If the requested Ritz pairs do not meet the true residual tolerance,
// eigenpairs are still returned but result.converged is false.
[[nodiscard]] SparseSymmetricEigenResult solve_sparse_symmetric_generalized_eigen(
    const CsrMatrix& stiffness,
    const CsrMatrix& mass,
    const SparseSymmetricEigenOptions& options = {});

} // namespace nova::linalg
