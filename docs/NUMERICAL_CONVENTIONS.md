# Numerical Conventions

## Coordinate and field conventions

- Internal calculations use SI units unless a module explicitly documents otherwise.
- Coordinates and field quantities are represented using IEEE 754 double precision in the initial CPU implementation.
- Arrays use zero-based indexing.
- Mesh and DOF identifiers are stable within a model revision; persistent project IDs are managed separately.
- A physics interface must document the units and sign convention of every field, source, flux, and boundary quantity.

## Weak-form convention

For the scalar diffusion equation

\[
-\nabla\cdot(k\nabla u)=f \quad \text{in }\Omega,
\]

with Dirichlet data on \(\Gamma_D\) and prescribed outward flux on \(\Gamma_N\), the weak form is

\[
\int_\Omega k\nabla u\cdot\nabla v\,d\Omega
=
\int_\Omega fv\,d\Omega
-
\int_{\Gamma_N} q_n v\,d\Gamma,
\]

where \(q_n = -k\nabla u\cdot n\) is defined as the outward physical flux. Individual physics modules must preserve this convention or explicitly declare a different one.

## Sparse matrix convention

- The initial sparse storage format is compressed sparse row (CSR).
- Column indices within each row are stored in ascending order.
- Duplicate assembly entries are summed.
- Non-finite matrix and vector values are rejected.
- Matrix dimensions and vector lengths are checked at API boundaries.

## Solver convention

The initial conjugate-gradient implementation is intended for symmetric positive-definite systems. It checks for non-positive curvature during iteration but does not constitute a formal proof that the input matrix is SPD.

Convergence is declared when

\[
\|r_k\|_2 \leq \max(\text{atol}, \text{rtol}\,\|r_0\|_2).
\]

A solver result must expose convergence status, iteration count, initial residual, final residual, and a diagnostic message. Reaching the iteration limit is not reported as successful convergence.

## Boundary conditions

The first Poisson benchmark applies strong Dirichlet constraints through row-and-column elimination with right-hand-side correction. More general constraint systems, periodic constraints, multipoint constraints, and Lagrange multipliers will be designed separately.

## Reproducibility

Every reported benchmark must record:
- NOVA revision.
- Compiler and build configuration.
- Mesh definition and element order.
- Material coefficients and units.
- Boundary and initial conditions.
- Solver type, tolerances, and iteration limit.
- Norm definitions and reference solution.
