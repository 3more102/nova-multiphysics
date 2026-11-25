# Sparse 3D PEC Maxwell Cavity Modes

## Problem

The lossless isotropic PEC cavity discretization solves the generalized symmetric
eigenproblem

\[
K e = \lambda M e, \qquad \lambda = k_0^2,
\]

where \(K\) is the Nédélec curl-curl matrix, \(M\) is the positive edge-mass
matrix, and \(f=c\sqrt{\lambda}/(2\pi)\). The curl-curl matrix is positive
semidefinite, not positive definite: discrete gradients are exact zero modes.
Applying a smallest-eigenvalue method without removing that kernel returns
static gradient fields instead of the first electromagnetic cavity modes.

## Exact discrete-gradient deflation

For the PEC boundary, nodal scalar potentials vanish on boundary nodes. Let
\(G\) map the interior H1 nodal values to oriented free-edge circulations.
For a globally oriented edge \((a,b)\), its gradient coefficient is
\(\phi_b-\phi_a\). The assembled edge-element complex satisfies

\[
K G = 0.
\]

The solver constructs \(G\) directly from mesh topology and the free-edge
ordering, then forms the sparse Gram matrix

\[
A_G = G^T M G.
\]

The M-orthogonal projection away from the exact-gradient space is

\[
P_M x = x-G A_G^{-1}G^T Mx.
\]

The projection system is solved with conjugate gradients and Jacobi
preconditioning. This removes the known discrete-gradient kernel without
building a dense nullspace basis. The number of interior nodal potentials is
reported as the exact-gradient nullity. Nontrivial topology can introduce
additional harmonic zero modes; zero Ritz values observed in the computed
subspace are counted separately and filtered from the positive spectrum.

## Positive shift-invert

For a positive shift \(\sigma\), define

\[
T_\sigma=(K+\sigma M)^{-1}M.
\]

For every generalized eigenpair, \(T_\sigma e=\mu e\), with

\[
\lambda=\frac{1}{\mu}-\sigma.
\]

When \(K\) is positive semidefinite and \(M\) is positive definite,
\(K+\sigma M\) is positive definite, so each operator application is an
SPD sparse CG solve. The largest positive transformed Ritz values correspond
to the smallest original eigenvalues. The exact-gradient complement is
invariant under this operator in exact arithmetic; the implementation
reprojects each transformed vector to control iterative-solve drift.

The implementation uses deterministic block starting vectors, two-pass
M-orthogonalization, a symmetric projected Rayleigh-Ritz problem, and true
residual checks in the original pencil. Only the projected matrix is dense,
with a hard dimension limit of 512. The full edge operators remain CSR.

## API and controls

Use `solve_maxwell_cavity_shift_invert_3d(mesh, matrices, options)`. Important
controls are:

| Option | Default | Meaning |
| --- | ---: | --- |
| `max_modes` | 6 | Number of lowest positive modes requested |
| `max_subspace` | 64 | Krylov/Ritz dimension, at most 512 |
| `shift_scale` | 0.05 | Multiplier for max(diag(K)/diag(M)) |
| `relative_tolerance` | 1e-8 | Required original-pencil residual |
| `shift_solve_relative_tolerance` | 1e-11 | Inner shifted-system CG tolerance |
| `projection_relative_tolerance` | 1e-11 | Gradient projection CG tolerance |
| `nullspace_relative_tolerance` | 1e-10 | Zero-mode cutoff relative to the diagonal spectral scale |

The result reports the shift, subspace size, operator applications, inner CG
iterations, exact-gradient nullity, additional zero Ritz modes detected,
frequencies, edge fields, and relative residuals. PEC edge coefficients are
returned as zero in the full global edge ordering.

## Verification

The automated tests compare the first three sparse modes against the existing
dense generalized-eigen reference on unit-cube meshes, check original-pencil
residuals, M-normalization and mutual M-orthogonality, verify gradient
orthogonality, check exact PEC zeros, and test the expected dielectric scaling
\(\lambda\propto\epsilon_r^{-1}\), \(f\propto\epsilon_r^{-1/2}\).

The dense solver is a small-matrix reference, not an independent experimental
truth source. The analytical PEC cube limit is \(k_0^2=2\pi^2\) for the
lowest family; mesh-convergence studies remain necessary.

## Current limitations

- This is an initial sparse shift-invert implementation, not a production
  Maxwell eigensolver qualification.
- Inner systems use Jacobi-preconditioned CG. Algebraic multigrid, robust
  incomplete-factorization fallback, parallel sparse kernels, and GPU execution
  are not implemented here.
- The projected Ritz problem is bounded to 512 dimensions. Large meshes may
  require restart/locking and a more scalable eigensolver backend.
- The exact-gradient construction handles the H1_0 image. Additional harmonic
  fields in multiply connected domains are detected only when represented in
  the current Krylov subspace; users must increase the subspace if the solver
  returns fewer than the requested positive modes.
- The automatically chosen shift is a scale-based starting point. Difficult
  meshes may require tuning the shift and inner CG tolerances.
- Passing these tests does not establish commercial readiness, general
  anisotropic-material support, lossy eigenmodes, PML support, or validated
  interoperability with external commercial solvers.
