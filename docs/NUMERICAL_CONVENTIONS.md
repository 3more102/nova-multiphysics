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

## Triangle quality and uniform mesh refinement

The 2D quality metric uses each triangle's positive area \(A\) and edge
lengths \(a,b,c\):

\[
q_\mathrm{mean}=\frac{4\sqrt{3}A}{a^2+b^2+c^2}.
\]

An equilateral triangle has \(q_\mathrm{mean}=1\); a degenerate shape
approaches zero. The smallest interior angle is computed using
\(\operatorname{atan2}(|u\times v|,u\cdot v)\) in radians.
Signed area preserves the source mesh orientation; total area sums
absolute areas.

Uniform 1-to-4 refinement creates one globally shared node per unique
original edge, then splits each triangle into three corner triangles
and one central triangle. The refinement is conforming on a valid
manifold triangular mesh; all source triangles retain orientation.
When exterior boundary edges are provided, each is split into two
segments with the same directed traversal. Original vertex values
are preserved and new midpoint temperatures for a P1 nodal field
are interpolated as the arithmetic mean of endpoint values.

This refinement is geometric and uniform, not a residual-based error
estimator or selective adaptive refinement. The reported quality is a
shape metric, not an independent solution-accuracy guarantee.

## H1 solution transfer and thermal remeshing

A conforming uniform refinement preserves all coarse vertices and adds
one mesh vertex at the midpoint of every unique coarse edge. For P1,
the new vertex value is the average of the two endpoint values. For P2,
the new vertex value is the existing coarse edge-midpoint DOF.

P2 transfer also evaluates the original quadratic finite-element
polynomial at every midpoint of a fine-mesh edge. When a fine edge lies
on a coarse-element interface, both parent-element traces are checked
for continuity before the value is accepted. This preserves the
piecewise-quadratic field exactly up to floating-point evaluation error.

The transient-heat remeshing helper applies the same conforming mesh
refinement, transfers an optional P1/P2 restart vector, and splits
Dirichlet, Neumann and Robin boundary edges into two directed segments.
Their original physical-coordinate callbacks are retained unchanged.
This is uniform remeshing; it does not yet implement error-driven local
refinement or transfer arbitrary user-defined physical groups.

## Local adaptive refinement and spatial error indicators

The one-level adaptive triangle refiner marks selected parent cells for
red refinement into four children. Every edge of a marked cell receives
one globally shared midpoint. Adjacent unmarked cells are green-refined
with one-, two-, or three-edge transition patterns, so neighboring cells
share the same split edges and no hanging nodes are introduced. Each
child retains its parent-cell index. A supplied boundary edge remains
one directed segment when untouched and becomes two directed segments
when bisected.

The spatial indicator uses area-weighted patch averaging of element
gradients to recover a continuous nodal gradient. It compares the P1
interpolation of that recovered gradient with the P1/P2 finite-element
gradient using the conductivity-weighted element energy norm:

eta_K^2 = integral_K k * |grad(T_recovered) - grad(T_h)|^2 dOmega.

Dörfler bulk marking selects the smallest deterministic set of cells,
ordered by descending indicator, whose squared indicators reach the
configured fraction of the global squared indicator sum. The recovery
indicator is a practical marking heuristic; it is not a guaranteed
upper bound on the true discretization error and does not replace
independent convergence studies.

The reference example demonstrates a complete single-cycle thermal
adaptation: solve on a coarse P2 mesh, recover the conductivity-weighted
gradient indicator, apply Dörfler marking, transfer the accepted P2 state
to the locally refined mesh, and continue the transient solve at the
coarse solution's final physical time.

Selective transient-heat remeshing transfers the optional P1/P2 restart
state, preserves spatial and time-dependent callbacks, and splits only
those Dirichlet, Neumann, or Robin edges whose parent edges were bisected.
This first adaptive refiner performs one conforming refinement level per
call; repeated marking cycles and more advanced grading/quality controls
remain application-level work.

## Tetrahedral 3D mesh and P1 H1 convention

The initial 3D geometry kernel uses four-node linear tetrahedra. The
reference tetrahedron has vertices (0,0,0), (1,0,0), (0,1,0), and
(0,0,1), with shape functions N0=1-r-s-t, N1=r, N2=s, N3=t.
The physical gradient is transformed by the inverse transpose of the
signed element Jacobian; the absolute determinant is used for volume
integration, while the signed determinant preserves orientation.

Tetrahedron validation rejects repeated/duplicate cells, invalid node
indices, non-finite coordinates, near-zero signed volume relative to
the cube of the longest edge, unused vertices, non-manifold faces and neighbors
occupying the same side of a shared face. Both element
orientations are accepted. Structured boxes use six tetrahedra per
hexahedral cell with a globally consistent face triangulation.

The tetrahedron mean-ratio quality is
q = 12 * (3V)^(2/3) / sum_i(edge_length_i^2),
where V is positive tetrahedron volume and the sum includes all six
edges. A regular tetrahedron has q=1. The minimum interior dihedral
angle is reported in radians. Edge vectors are divided by the longest edge
before forming unit face normals. The interior angle is acos(-n1 dot n2). Boundary faces are extracted once and
oriented outward using the opposite tetrahedron vertex.

The first 3D H1 space is continuous P1 with one global degree of freedom
per mesh vertex. The scalar diffusion kernel assembles
\(-\nabla\cdot(k\nabla u)=f\) on affine tetrahedra, using a positive
four-point degree-two tetrahedron rule. It supports spatially varying
positive isotropic diffusivity, volumetric source terms, strong Dirichlet
vertex/face values, outward Neumann fluxes, and Robin transfer on validated
exterior triangular faces. Three-point degree-two face quadrature integrates
linear Neumann loads and quadratic Robin matrices. Unassigned faces are
insulated. On each connected component without Dirichlet or positive Robin
anchoring, the integrated source must equal the integrated outward Neumann
flux within a relative tolerance of 1e-10 against the assembled absolute
load. The solver fixes the constant nullspace by imposing a temporary vertex
gauge. The physical field in `values` is shifted to zero volume-weighted mean on
that component. The original iterate in `linear_solver.solution` and its true
residual describe the same pinned system.
Anisotropic tensors remain unsupported.

The 3D diffusion kernel is P1 only. Higher-order tetrahedral spaces and 3D adaptive refinement remain open.

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

For a candidate interval of width $\Delta t$, the 2D adaptive driver takes one full step and two half steps from the same initial state. Use $p=1$ for Backward Euler and $p=2$ for Crank–Nicolson. Let $T_h$ be the endpoint after two half steps and $T_f$ the endpoint after one full step.

$$
e_i=\frac{T_{h,i}-T_{f,i}}{2^p-1},\qquad
r_i=\frac{e_i}{\mathrm{atol}+\mathrm{rtol}\max(|T_{n,i}|,|T_{h,i}|)}.
$$

The driver accepts when $\sqrt{\sum_i r_i^2/N}\leq 1$. It stores the result of the two half steps and combines their energy diagnostics. A rejected trial leaves the accepted time and field intact. The controller uses the local-error exponent $1/(p+1)$ and obeys the time interval, step bounds and attempt limit.

This is a local time-error indicator. Sharp transients and stiff dynamics can require tighter controls. The 2D controller on this branch reassembles its volume operators on each full or half trial.

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


## Transient P1 tetrahedral heat

For volume capacity $c=\rho c_p>0$ and isotropic conductivity $k>0$,
$M_{ij}=\int cN_iN_j$ and $K_{ij}=\int k\nabla N_i\cdot\nabla N_j$.
On an affine tetrahedron with constant $c$, $M^e_{ij}=cV(1+\delta_{ij})/20$.
The positive four-point degree-two rule integrates this mass exactly.
Variable coefficients are sampled at the same quadrature points. Their products with basis functions can exceed the quadrature degree.

On a Robin face let $R_{ij}=\int hN_iN_j$, with $h\ge0$.
Let $F_i=\int QN_i-\int_{\Gamma_N}q_nN_i+\int_{\Gamma_R}hT_{amb}N_i$.
The outward flux is $q_n=-k\nabla T\cdot n$. Each boundary face accepts one condition, independent of the vertex permutation. Unassigned faces are insulated.

For $L(t)=K+R(t)$ the shared time integrator uses

$$
[M+\theta\Delta t L_{n+1}]T_{n+1}
=M T_n+\Delta t[(1-\theta)(F_n-L_n T_n)+\theta F_{n+1}].
$$

Backward Euler ($\theta=1$) evaluates forcing only at the new endpoint. Crank–Nicolson ($\theta=1/2$) uses both endpoints. Initial/restart values are projected onto the prescribed start-time trace. At each subsequent endpoint the identity rows begin at their exact prescribed temperatures and receive zero CG updates. An insulated transient component retains its absolute temperature through the capacity term.

Heat content is $H=\mathbf{1}^{T}MT$. Source and boundary exchanges are integrated over each step, positive into the domain. Dirichlet energy is the sum of the constrained-row residuals of the unmodified heat equation. The reported imbalance is the change in $H$ minus all four signed energy inputs. The manufactured mixed-face test checks each exchange against an independent analytical surface or volume integral.
