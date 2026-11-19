# Conservative face-integrated heat reconciliation (experimental)

The steady P1 tetrahedral solver enforces the finite-element weak equations.
Raw cell-average conductive fluxes are generally **not** continuous across
material interfaces, nor conservative on each individual tetrahedron. The
separate global FEM energy report checks a global budget, not cell-by-cell flux
equilibrium.

The C++20 SDK function
reconcile_p1_face_heat_rates_3d(problem, temperatures) reconstructs **one
integrated heat rate in W per mesh face**. This is an optional postprocessor;
it does not modify the stiffness matrix, solved temperatures or original P1
gradient/flux diagnostics.

## Definition

For every tetrahedron K, an integrated outward rate H(K,F) is defined for
each of its four faces F. On an interior face shared by K and L,

    H(K,F) = -H(L,F).

The prescribed Neumann rate is the same degree-two triangular integral used
in the 3D thermal solver. Robin rate is the integral of h(T_h - T_ambient),
with T_h evaluated from the P1 face interpolation. Unassigned exterior faces
are insulated, i.e. zero. At an exterior Dirichlet face, the heat rate is an
unknown reaction.

The independent element target is the degree-two tetrahedral integral of
the sum of the spatial source and any per-tetrahedron source density.
Let H0_F be the raw P1 rate on an exterior Dirichlet face or the average of
the opposite outgoing raw rates from both adjacent tetrahedra on an interior
face. With positive face areas A_F, NOVA minimizes

    sum over adjustable F of ((H_F - H0_F) / A_F)^2

subject to the exact signed cell constraints

    sum over F on boundary(K) H(K,F) = integral_K Q dV.

The discrete optimality equations are a weighted dual-cell graph Laplacian.
They are assembled into NOVA's CSR matrix and solved by the existing CG
implementation. Pure-Neumann components have a one-dimensional multiplier
nullspace: NOVA checks net source/flux compatibility, pins only one dual
multiplier, and never introduces a fake physical face rate. This auxiliary
pin is unrelated to a Dirichlet temperature boundary.

## Interpreting the result

- faces contains one canonical (sorted-vertex) record per geometrical face.
  For an interior face, outward_heat is positive leaving first_cell and
  entering second_cell. For an exterior face, it is positive leaving the domain.
- generated_heat_per_cell and cell_balance_defect report the local budget in W.
- total_exterior_outward_heat and total_generated_heat should agree for an
  admissible balanced model, subject to floating-point solver tolerances.
- max_absolute_correction measures how much the raw P1 integrated rate had
  to change; correction_l2_flux_density is the norm of per-face corrections
  divided by area.

**Scope:** This is conservative *only in the face-integrated sense*. It is
not an H(div)-conforming finite-element vector field, a normal flux function
defined continuously at every point of each interface, a pointwise PDE
residual, a guaranteed discretization-error estimator, or a verified
physical constitutive reconstruction. Additional accuracy, convergence and
experimental validation are needed before engineering use.

Explicit node-only Dirichlet constraints are deliberately rejected because
their concentrated finite-element reactions cannot be assigned uniquely to
boundary faces. Give physical exterior Dirichlet face conditions instead.
The input field must agree with its prescribed temperatures. For unanchored
connected components, net volumetric generation and prescribed boundary
flux must be compatible.

## Reproduce

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel
    ctest --test-dir build --output-on-failure

The unit suite checks an exact two-material constant-flux patch, positive
volume generation, mixed Robin/Neumann conditions, orientation invariance,
heterogeneous SPD tensors with additive regional sources, floating component
compatibility and invalid data. This SDK endpoint has no standalone CLI
switch or VTK flux-field export yet.
