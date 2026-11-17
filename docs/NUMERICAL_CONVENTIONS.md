# Numerical Conventions

## Coordinate and field conventions

- Internal calculations use SI units unless a module explicitly documents otherwise.
- Coordinates and field quantities are represented using IEEE 754 double precision in the initial CPU implementation.
- Arrays use zero-based indexing.
- Mesh and DOF identifiers are stable within a model revision; persistent project IDs are managed separately.
- A physics interface must document the units and sign convention of every field, source, flux, and boundary quantity.

## Scalar diffusion weak form

For

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

where \(q_n=-k\nabla u\cdot n\) is the outward physical flux. The current 2D implementation supports volumetric source terms and nodal Dirichlet conditions; the Neumann term is documented for the future boundary-facet implementation.

## Reference triangle and P1 basis

The reference triangle has vertices \((r,s)=(0,0),(1,0),(0,1)\). Its linear shape functions are

\[
N_1=1-r-s,\qquad N_2=r,\qquad N_3=s.
\]

The current integration rule uses six positive symmetric points with reference weights approximately \(0.111690794839005\) (three points) and \(0.054975871827661\) (three points). The weights sum to the reference-triangle area \(1/2\) and the rule integrates polynomials through degree four exactly. This is sufficient for the current P1 scalar diffusion stiffness and for polynomial source terms through the corresponding assembled integrand degree.

Physical gradients are obtained from the element Jacobian. The signed determinant is retained in the gradient transformation; its absolute value is used for physical integration weights. Both triangle orientations are accepted, while degenerate and near-degenerate triangles are rejected.

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

The current 1D and 2D prototype paths apply strong Dirichlet constraints through row-and-column elimination with right-hand-side correction. General constraint systems, periodic constraints, multipoint constraints, Neumann and Robin boundary integration, and Lagrange multipliers remain future work.

## Reproducibility

Every reported benchmark must record:
- NOVA revision.
- Compiler and build configuration.
- Mesh definition and element order.
- Material coefficients and units.
- Boundary and initial conditions.
- Solver type, tolerances, and iteration limit.
- Norm definitions and reference solution.
