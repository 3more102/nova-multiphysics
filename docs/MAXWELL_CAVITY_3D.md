# Lossless PEC cavity eigenmodes (bounded reference solver)

The existing orientation-consistent lowest-order edge basis is used for
K e = lambda M e, where K_ij = integral(mu_r^-1 curl N_i dot curl N_j),
M_ij = integral(epsilon_r N_i dot N_j), lambda = k0^2. All exterior
tangential electric DOFs vanish. Frequency = c sqrt(lambda)/(2 pi).
Real positive scalar cell materials only; dispersive, complex or negative
permittivity is not supported by this eigenproblem.

The general symmetric solver accepts a real symmetric sparse pair, forms
M=L L^T by Cholesky, and solves the symmetric dense L^-1 K L^-T problem
with cyclic Jacobi rotations. Eigenvectors are transformed by L^-T and
normalized in the M inner product. A scaled residual is reported for every
mode. Convergence is checked against the off-diagonal Frobenius norm and
residuals. This is limited to 512 unknowns (O(n^3) time/O(n^2) storage);
it is not a scalable shift-invert eigensolver.

The discrete gradient kernel is filtered by an explicit relative
spectral threshold (default 1e-10 of the maximum absolute eigenvalue).
Negative eigenvalues beyond this cutoff are rejected. Nullity is reported.
This numerical filter is not a general topology-aware spurious-mode proof
and may remove genuine very small positive eigenvalues in ill-scaled
problems. Degenerate modes are returned as a basis, not individually tracked.

## Independent reference and verification

For a unit PEC vacuum cube the first three physical eigenvalues are
2*pi^2 (19.7392088 m^-2). The gradient nullspace has (n-1)^3 interior
nodal potentials for this structured connected cube. Independent SciPy
generalized dense eigensolutions with numerical mass integration:

| Cells | Free edges | Nullity | First lambda | Second/third lambda | First frequency relative error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 48 | 26 | 1 | 17.06363423 | 19.64300762 | -7.024% |
| 162 | 117 | 8 | 18.43103487 | 20.02505572 | -3.370% |
| 384 | 316 | 27 | 18.96183604 | 19.94375703 | -1.989% |

These are measured independent reference results, not C++ results until
CI executes. The automated C++ tests compare all three eigenvalues against
the reference to 2e-7 absolute, check the nullity, scaled residual <1e-10,
zero exterior edge coefficients, and decreasing first-frequency error.
A separate check verifies that quadrupling epsilon halves the frequency.
The generic solver tests a known generalized spectrum and M-orthogonality
and rejects nonsymmetric or indefinite-mass inputs.

Reproduce using the existing CMake/CTest workflow or
python verification/verify_maxwell_cavity_reference.py (NumPy/SciPy).
No new C++ dependency; optional Python verification only.
