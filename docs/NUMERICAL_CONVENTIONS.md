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

Neumann conditions supply outward physical flux \(q_n\) and contribute \(-\int_{\Gamma_N}q_n N_i\,d\Gamma\) to the right-hand side. Robin conditions contribute \(\int_{\Gamma_R}h u v\,d\Gamma\) to the left-hand side and \(\int_{\Gamma_R}h u_\infty v\,d\Gamma\) to the right-hand side. Three-point Gauss-Legendre quadrature is used on P1 and P2 boundary edges (exact for polynomials through degree five). Neumann and Robin conditions cannot be assigned to the same edge.

## Transient heat diffusion

The transient scalar heat equation is

\[
\rho c\frac{\partial T}{\partial t}-\nabla\cdot(k\nabla T)=Q
\]

with spatially varying, time-independent, strictly positive volumetric heat capacity \(\rho c\) and isotropic conductivity \(k\). For H1 basis \(N_i\),

\[
M_{ij}=\int_\Omega\rho cN_iN_j\,d\Omega,\quad
K_{ij}=\int_\Omega k\nabla N_i\cdot\nabla N_j\,d\Omega,\quad
b_{Q,i}(t)=\int_\Omega Q(t)N_i\,d\Omega.
\]

On an exterior edge with outward conductive heat flux \(q_n=-k\nabla T\cdot n\), the Neumann load is \(b_{N,i}=-\int_{\Gamma_N}q_nN_i\,d\Gamma\). Positive outward heat flux is an energy *loss*. For Robin convection \(q_n=h(T-T_\infty)\), the weak form adds \(R_{ij}(t)=\int_{\Gamma_R}hN_iN_j\,d\Gamma\) to the left and \(b_{R,i}(t)=\int_{\Gamma_R}hT_\infty N_i\,d\Gamma\) to the right, with \(h\ge 0\). Overlapping boundary conditions on the same edge are rejected. Unassigned edges are insulated.

The mass and stiffness operators use positive degree-four six-point reference-triangle quadrature. Natural boundary contributions use three-point Gauss–Legendre edge quadrature. P2 edge-midpoint DOFs are shared globally and included in natural boundary integrals.

The theta time integrator chooses \(\theta=1\) for Backward Euler and \(\theta=1/2\) for Crank–Nicolson:

\[
[M+\theta\Delta t(K+R_{n+1})]T^{n+1} =
[M-(1-\theta)\Delta t(K+R_n)]T^n+
\Delta t[\theta F_{n+1}+(1-\theta)F_n],
\]

where \(F=b_Q+b_N+b_R\). Strong prescribed temperatures are imposed at each new time level using symmetric elimination with a corrected right-hand side. Initial essential boundary DOFs are projected onto their prescribed temperatures at the configured start time.

The per-step *integrated* discrete heat budget is

\[
E_{n+1}-E_n =
E_Q+E_N+E_R+E_D+\varepsilon,\qquad
E_n=\int_\Omega\rho cT_h^n\,d\Omega.
\]

The energy terms integrate their respective weighted theta-method source and boundary contributions over one time step. The Dirichlet term \(E_D\) is calculated by summing constrained finite-element residual entries; it is a *discrete reaction*, not an independently reconstructed normal flux. The energy-balance residual \(\varepsilon\) checks an algebraic identity and iterative solver error, **not** independent physical accuracy. Independent numerical verification is still needed.

Backward Euler has first-order temporal accuracy; Crank–Nicolson has second-order temporal accuracy for sufficiently smooth solutions. Both are implicit, but Crank–Nicolson is not L-stable and large time steps may be inaccurate or nonmonotone. The fixed-step API supports explicit physical start times and an exact initial H1 DOF restart vector, including P2 edge-midpoint values. Time-dependent material coefficients remain unsupported.

### Adaptive time-stepping convention

For a candidate interval of width (Delta t), the adaptive driver performs one full theta-method step and two half steps from the **same initial H1 state**. Define (p=1) for Backward Euler or (p=2) for Crank–Nicolson. Its local error indicator at each DOF is

[
e_i=rac{T_{2	imesmathrm{half},i}-T_{mathrm{full},i}}{2^p-1},
qquad
r_i=rac{e_i}{mathrm{atol}+mathrm{rtol}max(|T_{n,i}|,|T_{2	imesmathrm{half},i}|)}.
]

It accepts a trial when (sqrt{sum_i r_i^2/N}le 1). On acceptance it stores the **two-half-step** solution, not Richardson extrapolation, and aggregates their energy diagnostics. On rejection the physical state does not advance. A bounded controller adjusts the next step using the local-error exponent (1/(p+1)), while enforcing the configured step bounds, time interval and attempt limit.

This is a local DOF-RMS temporal error criterion, not a global error bound or a spatial error estimate. Very sharp transients, stiff dynamics, and discontinuous coefficients can require tighter controls. The current reference implementation reassembles its matrices on each full/half trial and is intended for small benchmarks rather than high-performance production use.

## Nonlinear transient heat conduction

The nonlinear P1/P2 H1 prototype solves

\[
\rho c(x,y)\frac{\partial T}{\partial t}
-\nabla\cdot(k(x,y,T)\nabla T)=Q(x,y,t).
\]

For Backward Euler and a test function \(N_i\), the unconstrained discrete residual at time \(t_{n+1}\) is

\[
R_i(T^{n+1})=
\int_\Omega
\left[
\frac{\rho c}{\Delta t}(T_h^{n+1}-T_h^n)N_i
+k(x,y,T_h^{n+1})\nabla T_h^{n+1}\cdot\nabla N_i
-Q(t_{n+1})N_i
\right]d\Omega
+\int_{\Gamma_N}q_nN_i\,d\Gamma
+\int_{\Gamma_R}h(T_h^{n+1}-T_\infty)N_i\,d\Gamma.
\]

The *consistent* tangent is

\[
\frac{\partial R_i}{\partial T_j}=
\int_\Omega \left[
\frac{\rho c}{\Delta t}N_iN_j
+k\nabla N_j\cdot\nabla N_i
+\frac{\partial k}{\partial T}N_j
  \nabla T_h^{n+1}\cdot\nabla N_i
\right]d\Omega
+\int_{\Gamma_R}hN_iN_j\,d\Gamma.
\]

Here \(k>0\), \(\rho c>0\), and the exact constitutive derivative \(\partial k/\partial T\) must be supplied. The Robin coefficient \(h\) may depend on position and time but **not** on temperature in this version. Strong Dirichlet constraints overwrite residual rows with \(T_i-g_i\) and tangent rows with the corresponding identity row. Unconstrained columns are retained in free equations; this is the exact nonsymmetric Jacobian for the complete unknown vector.

Damped Newton iteration uses the existing Armijo line search and either GMRES or BiCGSTAB, not Conjugate Gradient. Its state is advanced only after a step converges. Nonconvergence is reported as a failure rather than a valid transient solution. This solver uses fixed Backward Euler steps only; the linear adaptive Crank–Nicolson controller does not currently apply to temperature-dependent material laws.

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
- **Jacobi:** diagonal preconditioning \(M^{-1}=\operatorname{diag}(A)^{-1}\). CG requires every diagonal entry to be finite and strictly positive. BiCGSTAB and GMRES require every diagonal entry to be finite and nonzero.

GMRES uses right Jacobi preconditioning, modified Gram-Schmidt Arnoldi orthogonalization with a second reorthogonalization pass, and Givens rotations. Its restart length limits the Krylov subspace memory.

Before returning success, the iterative linear solvers recompute the true residual \(b-Ax\). A recursive or projected residual that meets tolerance is not sufficient if the recomputed residual does not.

## Nonlinear Newton convention

The Newton kernel solves \(F(x)=0\) using the Jacobian \(J(x)=dF/dx\) and a linear correction

\[
J(x_k)\,\Delta x_k=-F(x_k).
\]

A backtracking line search tests \(x_{k+1}=x_k+\lambda\Delta x_k\) and accepts a finite trial when

\[
\|F(x_k+\lambda\Delta x_k)\|_2
\leq
(1-c\lambda)\|F(x_k)\|_2,
\]

where \(0<c<1\) is the Armijo constant and \(\lambda\) is reduced by one half after a rejected trial. The nonlinear iteration is successful only when the residual norm reaches the configured tolerance. A small step without residual convergence is reported as stagnation, not success.

## Nonlinear 1D diffusion verification

The current nonlinear FEM benchmark uses the constitutive coefficient

\[
k(u)=1+u^2,\qquad k'(u)=2u.
\]

For a P1 field \(u_h\), the element residual is assembled from

\[
R_i^e=\int_{\Omega_e} k(u_h)u_h' N_i'\,dx
-\int_{\Omega_e}fN_i\,dx,
\]

and the consistent tangent is

\[
J_{ij}^e=\int_{\Omega_e}
\left[k'(u_h)N_j u_h' N_i'
+k(u_h)N_j'N_i'\right]dx.
\]

Three-point Gauss-Legendre quadrature is used on each line element. The manufactured solution \(u(x)=\sin(\pi x)\) defines the source through \(f=-d[k(u)u']/dx\). Endpoint Dirichlet values are imposed strongly.

The nonlinear FEM benchmark is a verification example; it is not yet a general nonlinear constitutive-material API.

## Reproducibility

Every benchmark must record the NOVA revision, compiler, build configuration, mesh, element order, material coefficients, units, boundary conditions, solver type, preconditioner, tolerances, GMRES restart length, nonlinear line-search parameters, norm definitions, and reference solution.
