# Sparse generalized symmetric eigen solver (experimental)

## Purpose

`solve_sparse_symmetric_generalized_eigen` computes the largest algebraic Ritz values of

`K x = lambda M x`

without converting the full matrices to dense storage. It is intended for bounded reference problems and as a scalable foundation for modal FEM workflows.

## Numerical method

- Deterministic block starting vectors reduce the risk of missing repeated or clustered dominant eigenvalues.
- The generalized operator action is evaluated as `M^{-1} K q`, with each mass system solved using the existing Jacobi-preconditioned conjugate-gradient implementation.
- Two-pass modified Gram-Schmidt maintains orthogonality in the `M` inner product.
- When the Krylov queue reaches an invariant subspace before the projection limit, deterministic coordinate seeds explore the remaining M-orthogonal complement.
- A symmetric Rayleigh-Ritz projection is assembled on the generated subspace.
- The projected dense eigensystem is solved with NOVA's existing symmetric generalized eigensolver.
- Ritz vectors are reconstructed in the original sparse space, normalized with `x^T M x = 1`, and checked using the true sparse residual.

The solver returns the largest available algebraic Ritz pairs in descending order. A failed residual tolerance or insufficient returned pair count sets `converged=false`. A small residual certifies an approximate eigenpair, not its global spectral ordering: a bounded projection smaller than the full dimension can still miss an eigenvalue. The complement restart prevents premature termination when an invariant subspace exhausts the Krylov queue and unused projection capacity remains.

## Limits

- `K` must be real symmetric; `M` must be symmetric positive definite.
- Mass solves use Jacobi-preconditioned CG. Poorly conditioned or strongly graded mass matrices may require more iterations or a future stronger preconditioner.
- The projected dense Rayleigh-Ritz problem is limited to 512 dimensions.
- The current implementation uses a fixed deterministic starting block and does not yet implement thick restart, shift-invert, spectrum slicing, or a sparse direct factorization.
- A failed mass solve raises an exception. The algorithm does not silently continue with an inaccurate inverse action.
- This is an experimental largest-algebraic-eigenvalue solver, not a certified eigensolver for engineering release.

## Verification

The regression suite covers:
- A generalized diagonal pencil with a repeated dominant eigenvalue.
- A two-dimensional pencil whose deterministic initial vector lies exactly in a nondominant invariant subspace.
- Agreement with the dense generalized reference on a coupled SPD pencil.
- Explicit reporting of unconverged Ritz pairs.
- A 600-DOF sparse case, above the dense reference limit.
- Invalid dimensions, nonsymmetric operators, non-positive mass diagonals, and invalid options.

Build and run:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```
