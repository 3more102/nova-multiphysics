# Numerical Conventions

## Coordinate and field conventions
- Internal calculations use SI units unless a module explicitly documents otherwise.
- Coordinates and field quantities use IEEE 754 double precision in the initial CPU implementation.
- Arrays use zero-based indexing.
- A physics interface documents units and sign conventions for fields, sources, fluxes, and boundary quantities.
- Initial 2D triangular elements are defined in the XY plane.

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
\int_{\Gamma_N}q_n v\,d\Gamma,
\]

where \(q_n=-k\nabla u\cdot n\) is the outward physical flux. The current 2D implementation supports volumetric source terms and nodal Dirichlet conditions; natural boundary integration is future work.

## Reference triangle and P1 basis

The reference triangle has vertices \((r,s)=(0,0),(1,0),(0,1)\). Its linear shape functions are

\[
N_1=1-r-s,\qquad N_2=r,\qquad N_3=s.
\]

The assembly rule uses six positive symmetric points and integrates polynomials through degree four exactly.

Physical gradients are obtained from the element Jacobian. The signed determinant is retained in the gradient transformation; its absolute value is used for physical integration weights. Both triangle orientations are accepted, while degenerate and near-degenerate triangles are rejected using a scale-aware test.

## Error-norm quadrature

The L2 and H1 error calculator uses a 5×5 Gauss-Legendre product rule mapped to the reference triangle by a Duffy transformation. The resulting rule integrates polynomial integrands through total degree eight exactly, which is sufficient for the squared P1 interpolation error of the current quartic manufactured solution.

For numerical field \(u_h\) and exact field \(u\):

\[
\|u-u_h\|_{L^2(\Omega)}
=
\left(\int_\Omega (u-u_h)^2\,d\Omega\right)^{1/2},
\]

\[
|u-u_h|_{H^1(\Omega)}
=
\left(\int_\Omega |\nabla u-\nabla u_h|^2\,d\Omega\right)^{1/2},
\]

\[
\|u-u_h\|_{H^1(\Omega)}
=
\left(\|u-u_h\|_{L^2(\Omega)}^2+
|u-u_h|_{H^1(\Omega)}^2\right)^{1/2}.
\]

## Sparse matrix convention
- Initial sparse storage is compressed sparse row (CSR).
- Column indices within each row are stored in ascending order.
- Duplicate assembly entries are summed.
- Non-finite matrix and vector values are rejected.
- Matrix dimensions and vector lengths are checked at API boundaries.

## Solver convention

Conjugate Gradient is intended for symmetric positive-definite systems. It detects non-positive search-direction curvature but cannot prove that an arbitrary input matrix is SPD.

Supported preconditioners:
- **None:** standard unpreconditioned CG.
- **Jacobi:** diagonal preconditioning \(M^{-1}=\operatorname{diag}(A)^{-1}\). Every diagonal entry must be finite and strictly positive.

Convergence is defined by

\[
\|r_k\|_2 \leq \max(\text{atol},\text{rtol}\,\|r_0\|_2).
\]

Before returning success, the solver recomputes the true residual \(b-Ax\). If the recursive residual meets tolerance but the recomputed residual does not, the solver restarts from the true residual.

A solver result exposes convergence status, iteration count, initial residual, final true residual, and a diagnostic message. Reaching the iteration limit is not successful convergence.

## Boundary conditions

The current 1D and 2D prototype paths apply strong Dirichlet constraints through row-and-column elimination with right-hand-side correction. General constraint systems, periodic constraints, multipoint constraints, Neumann and Robin boundary integration, and Lagrange multipliers remain future work.

## Reproducibility

Every benchmark must record the NOVA revision, compiler, build configuration, mesh, element order, material coefficients, units, boundary conditions, solver type, preconditioner, tolerances, norm definitions, and reference solution.
