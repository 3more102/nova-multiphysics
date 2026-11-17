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

with Dirichlet data on \(\Gamma_D\) and prescribed outward physical flux \(q_n=-k\nabla u\cdot n\) on \(\Gamma_N\), the weak form is

\[
\int_\Omega k\nabla u\cdot\nabla v\,d\Omega
=
\int_\Omega fv\,d\Omega
-
\int_{\Gamma_N}q_n v\,d\Gamma.
\]

For a Robin condition, the outward flux is \(q_n=h(u-u_\infty)\), with non-negative transfer coefficient \(h\). The weak form adds \(\int_{\Gamma_R}h u v\,d\Gamma\) to the left-hand side and \(\int_{\Gamma_R}h u_\infty v\,d\Gamma\) to the right-hand side.

## Reference triangle and P1 basis

The reference triangle has vertices \((r,s)=(0,0),(1,0),(0,1)\). Its linear shape functions are

\[
N_1=1-r-s,\qquad N_2=r,\qquad N_3=s.
\]

The assembly rule uses six positive symmetric points and integrates polynomials through degree four exactly. Physical gradients are obtained from the element Jacobian. The signed determinant is retained in the gradient transformation; its absolute value is used for physical integration weights. Both triangle orientations are accepted, while degenerate and near-degenerate triangles are rejected using a scale-aware test.

## Boundary-edge integration

Neumann conditions supply outward physical flux \(q_n\) and contribute \(-\int_{\Gamma_N}q_n N_i\,d\Gamma\) to the right-hand side. Robin conditions contribute \(\int_{\Gamma_R}h u v\,d\Gamma\) to the left-hand side and \(\int_{\Gamma_R}h u_\infty v\,d\Gamma\) to the right-hand side. Two-point Gauss quadrature is used on linear edges. Neumann and Robin conditions cannot be assigned to the same edge.

## Error-norm quadrature

The L2 and H1 error calculator uses a 5×5 Gauss-Legendre product rule mapped to the reference triangle by a Duffy transformation. The resulting rule integrates polynomial integrands through total degree eight exactly, which is sufficient for the squared P1 interpolation error of the current quartic manufactured solution.

## Sparse matrix convention
- Initial sparse storage is compressed sparse row (CSR).
- Column indices within each row are stored in ascending order.
- Duplicate assembly entries are summed.
- Non-finite matrix and vector values are rejected.
- Matrix dimensions and vector lengths are checked at API boundaries.

## Solver convention

Conjugate Gradient is intended for symmetric positive-definite systems. BiCGSTAB and restarted GMRES are intended for general nonsymmetric square systems. All methods report breakdown conditions; none proves that an arbitrary matrix is suitable for the selected algorithm.

Supported preconditioners:
- **None:** no preconditioning.
- **Jacobi:** diagonal preconditioning \(M^{-1}=\operatorname{diag}(A)^{-1}\). CG requires every diagonal entry to be finite and strictly positive. BiCGSTAB and GMRES require finite nonzero entries.

GMRES uses right Jacobi preconditioning, modified Gram-Schmidt Arnoldi orthogonalization with a second reorthogonalization pass, and Givens rotations. Its restart length limits Krylov-subspace memory.

Before returning success, the iterative solvers recompute the true residual \(b-Ax\). A recursive or projected residual that meets tolerance is not sufficient if the recomputed residual does not.

## Reproducibility

Every benchmark must record the NOVA revision, compiler, build configuration, mesh, element order, material coefficients, units, boundary conditions, solver type, preconditioner, tolerances, GMRES restart length, norm definitions, and reference solution.
