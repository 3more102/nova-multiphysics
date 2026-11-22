# 3D Electrostatics on Tetrahedral Meshes

## Model

The experimental P1 interface solves the electrostatic Poisson problem

\[
-\nabla\cdot(\epsilon\nabla V)=\rho_{free},
\qquad E=-\nabla V,
\qquad D=\epsilon E.
\]

SI units are used: potential V [volt], permittivity epsilon [F/m], free charge density rho [C/m^3], electric field E [V/m], and electric displacement D [C/m^2].

The weak form reuses NOVA's validated scalar/tensor P1 tetrahedral diffusion assembly. Fixed node or exterior-face potentials are Dirichlet data. Prescribed exterior displacement is D dot n, with the same outward-flux sign convention as the diffusion solver. Unassigned exterior faces have zero normal displacement. Spatial charge and optional per-cell charge are additive; an omitted spatial charge callback means zero.

## Returned quantities

- Nodal electric potential.
- Per-cell constant electric field E = -grad(V).
- Per-cell quadrature-averaged displacement D = epsilon E.
- Stored field energy U = 1/2 integral(E dot D) dV.

The energy is evaluated in long double and checked for non-finite or materially negative density. Positive-definite permittivity implies E dot D >= 0.

## Fixed-electrode charge diagnostic (experimental)

Use `nova/physics/electrode_charge_3d.hpp` after a converged solution to
integrate the adjacent P1 tetrahedral electric displacement over selected
fixed-potential exterior triangles:

```cpp
auto charges = nova::physics::evaluate_electrode_charges_3d(
    problem, solved, electrode_faces, electrode_ids
);
```

`electrode_faces` and `electrode_ids` must have equal lengths. Each
nonnegative ID identifies one equipotential conductor boundary. All selected
vertices must be constrained by `potential_constraints` or `potential_faces`.
The routine rejects unconverged solutions, non-exterior faces, duplicates,
unfixed nodes, nonuniform group potentials and non-finite fields. It never
runs a second field solve.

For the dielectric domain's outward unit normal `n`, the integrated
displacement `integral(D dot n) dA` has units of coulombs. The charge on
the adjacent conductor is its **negative**. Each result reports both signs,
the electrode's potential, area and number of triangles. The reported groups
are sorted by electrode ID; input winding does not affect orientation.

In the included parallel-plate test (1 m² plates, 1 m separation, epsilon
= 2 F/m, V_low = 0 V, V_high = 10 V), the grounded electrode carries
-20 C and the high-potential electrode +20 C. Under those *charge-free*
two-electrode conditions, `C = Q_high/(V_high - V_low) = 2 F`, matching
`U = 100 J`.

**Accuracy warning:** This is a direct integration of *element-averaged*
`D`, not an equilibrated boundary flux, discrete electrode reaction, or
conservative surface-charge reconstruction. Coarse meshes, material jumps
or volume free charge can spoil charge balance. Capacitance from `Q/dV`
requires an appropriate charge-free two-terminal problem and independent
mesh-convergence verification.

## Verification

The regression suite includes:

1. A 1 m cube parallel-plate capacitor with epsilon = 2 F/m and 10 V separation: C = 2 F and U = 100 J.
2. An anisotropic SPD dielectric tensor with an affine exact potential, exact E and D, and U = 36.5 J.
3. A one-free-DOF tetrahedron with spatial plus per-cell charge density, giving V = 1 V and U = 1/12 J.
4. The same one-free-DOF tetrahedron with cell-only charge, verifying that an omitted spatial source is treated as zero.
5. A pure-Neumann affine field driven by prescribed D dot n, verifying outward-displacement sign, zero-volume-mean gauge and stored energy.
6. Rejection of negative permittivity, conflicting material representations and non-finite charge.

The installed SDK consumer also solves the capacitor benchmark and verifies the fixed-electrode charge API, checking the public headers and exported link symbols outside the source tree.

## Discrete weak-form electrode reaction (experimental)

The `evaluate_discrete_electrode_reactions_3d` public API is a **separate**
charge diagnostic for solved, fixed-potential conductor boundaries. Unlike
`evaluate_electrode_charges_3d`, which integrates the element-averaged
displacement over exterior triangular faces, this method uses the same
Galerkin volume weak form as the P1 solver:

\[
Q_g = \sum_{a\in E_g} (K V)_a
    = \sum_T \int_T \epsilon\nabla w_g\cdot\nabla V\,dT.
\]

Here `E_g` contains the boundary nodes belonging to electrode `g`, and
`w_g` is their P1 nodal indicator (one at those nodes, zero elsewhere).
For source-free, insulated models the result is `(K V)_g`. For models with
impressed charge and prescribed outward displacement it is instead the
**full Galerkin reaction**:

```text
Q_g = sum_{a in E_g} (K V - F)_a
    = integral epsilon grad(w_g) dot grad(V) dV
      - integral rho_free w_g dV
      + integral_N (D dot outward_normal) w_g dA
```

The volume and exterior-face loads use the same positive degree-two tetrahedron
and triangle quadratures as the stationary P1 FEM solver. Spatial and per-cell
charge are additive. The conductivity/permittivity term still uses
`-volume * grad(N_a) dot averaged_D`. All material representations already
supported by the solver remain available.

**Required conditions:** each selected electrode group must have its complete
exterior triangular contact supplied; no vertex may belong to two electrodes;
all selected vertices must have fixed potentials and one common equipotential
per group. Source and prescribed outward-displacement callbacks must be
finite and the per-cell charge array must match the tetrahedron count.
Neumann loads carry a positive sign in `Q=K V-F` because the FEM right-hand
side contains their *negative* outward flux.

The regression suite checks an independent Gauss-law balance for a constant
volume charge in a parallel-plate slab and for a prescribed nonzero side-face
displacement. On these manufactured geometries the sum of electrode reactions
matches the volume charge and prescribed boundary displacement consistently.
Such discrete consistency is not a certified continuum charge error.

For a linear, source-free, no-prescribed-flux system with insulated remaining
exterior, the discrete FEM identity `sum_g (V_g Q_g)=2*U` holds when all
nonzero Dirichlet electrodes are included. On the parallel-plate benchmark,
the two reaction charges are exactly approximately -20 C and +20 C (subject to
linear solver tolerance), and `V_high*Q_high=200 J=2*U`. A separate
variable-permittivity regression checks this energy identity on a coarse mesh.

**Caveat:** these are *discrete FEM reactions*, not exact physical surface
charges. They depend on the Galerkin approximation and mesh quality. This API
does not change the capacitance extraction matrix, solve floating-electrode
constraints, or certify physical charge conservation. It provides an
independent discretization-consistency observable for future work.

## Multi-terminal capacitance matrix (experimental)

Use `extract_capacitance_matrix_3d` with one group of exterior triangles per
conductor. It performs one electrostatic solve per electrode, setting that
electrode to 1 V and every other electrode to 0 V. Existing non-electrode
Dirichlet constraints and faces are treated as grounded references; existing
Dirichlet faces that exactly match an electrode group are replaced by the
basis voltage. The extraction excludes volume free charge and prescribed
displacement loads.

The returned electrode IDs are sorted. `matrix_farads` is row-major and
`matrix_farads[i*n+j]` is the conductor charge on electrode `i` for a
unit-voltage excitation on electrode `j`. The reported absolute and relative
reciprocity defects measure the asymmetry between `C_ij` and `C_ji`; the
matrix is not silently symmetrized. A single electrode requires a separate
grounded potential reference.

This first extraction path uses the same element-averaged displacement
boundary integral as the fixed-electrode diagnostic. It is therefore an
experimental estimate, not a discrete FEM electrode reaction or a locally
equilibrated charge. Check reciprocity and mesh convergence before relying on
the result. A separate terminal-algebra utility for prescribed floating charges is documented
below; a fully coupled floating-conductor FEM field solve remains future work.

The parallel-plate regression returns the Maxwell matrix
`[[2, -2], [-2, 2]] F` for `epsilon = 2 F/m`, unit area and unit gap.

## Weak-form discrete capacitance matrix (experimental)

`extract_discrete_capacitance_matrix_3d` uses the **same validated
unit-voltage basis FEM solves** as `extract_capacitance_matrix_3d`, but its
matrix entries are the electrode weak-form reactions `(K V)_electrode`
rather than raw outward surface fluxes. This implementation shares the basis
solver and matrix construction, so it does not run duplicate solves simply
to compare two charge-estimator methods.

The returned `CapacitanceMatrix3D` has the same electrode-ID ordering,
unsymmetrized row-major matrix, and reciprocity diagnostics. For a source-free
symmetric SPD dielectric and converged Galerkin solves, the discrete matrix
should be symmetric (within solver and floating-point tolerances). A two-plate
test checks `[[2,-2],[-2,2]] F`; a spatially varying-permittivity test checks
reciprocity and constant-voltage nullspace behavior.

Do not assume this matrix is the continuum Maxwell capacitance. It is based
on P1 Galerkin reaction charge and depends on mesh discretization. The
existing boundary-flux matrix remains available independently; comparing
both and repeating after mesh refinement helps identify flux-extrapolation
error. No hidden matrix symmetrization, voltage gauge, or physical calibration
is applied. Floating-terminal voltage inference can consume either matrix,
but a true nonlinear or charge-loaded conductor solve is not implemented.

## Floating-electrode terminal potentials (experimental)

After extracting a capacitance matrix, the *terminal-only* inverse
`solve_floating_electrode_potentials_3d` in
`nova/physics/floating_electrodes_3d.hpp` solves

```text
C_FF * V_F = Q_F - C_FK * V_K
```

where F indexes floating conductors with prescribed total charge Q, and K
indexes conductors with known voltage V. Every electrode ID must be assigned
exactly once. ID order in the input collections is arbitrary; results follow
the capacitance matrix's ascending electrode IDs. The API returns all terminal
voltages, estimated terminal charges from the original raw C, and the maximum
absolute floating-charge residual. It accepts an all-floating configuration
only when the capacitance block is nonsingular (for example, a separate
grounded reference exists). An isolated capacitance matrix with no reference
has a gauge nullspace, which the solver rejects rather than setting a hidden
ground. Singular or severely ill-conditioned blocks are also rejected.

For the two-terminal parallel plate benchmark with C = [[2,-2],[-2,2]] F,
fixing the first electrode to 0 V and prescribing +20 C on the second
returns 10 V on the second and -20 C on the first.

**Limitations:** This operation is a linear algebraic terminal inference, not
a new FEM field solve, not an independent electrostatic solution with floating
conductor boundary constraints, and not a model of charge in the dielectric.
As with the extracted matrix, source-free linear dielectrics and a valid
voltage reference are assumed. No raw-flux symmetrization or inferred gauge
is applied. Verify matrix reciprocity and mesh convergence, and independently
re-solve the FEM potential field using the computed terminal voltages if
spatial field quantities are required.

## Terminal-charge-to-field reconstruction (experimental)

The public `reconstruct_floating_electrostatic_field_3d` helper in
`nova/physics/floating_electrodes_3d.hpp` now completes a **source-free,
linear** terminal-to-field workflow:

1. Extract the unsymmetrized Maxwell capacitance matrix from electrode faces.
2. Infer the voltages of prescribed-charge floating terminals.
3. Reapply those electrode voltages as constant Dirichlet face data and solve
   the P1 electrostatic field once more.
4. Independently integrate the reconstructed electrode displacement fluxes and
   report the maximum difference against the requested floating charges.

For a unit parallel-plate example with `epsilon=2 F/m`, the low electrode
fixed to 0 V and the high electrode prescribed +20 C, the recovered high
terminal voltage is 10 V, the field is `E=(-10,0,0) V/m`, and the stored
energy is 100 J. Tests also cover reversed electrode-face winding, a single
floating terminal with an independent grounded exterior reference, and
rejection of source terms.

**Important conditions:** The input must have no volume free-charge callback,
no cell free-charge values, and no prescribed displacement-flux faces, because
the capacitance matrix excludes all such loads. Non-electrode original
Dirichlet faces/nodes are explicitly set to 0 V in the new solve, matching the
capacitance basis problems; existing electrode-face voltages are replaced by
the computed terminal values. The caller must supply a nonsingular capacitance
reference and keep electrode groups disjoint. Every reconstructed charge is
calculated from the same non-equilibrated P1 element displacement used for
capacitance extraction; a small residual is an **internal consistency check**,
not independent physical validation. Mesh convergence and reciprocity remain
caller responsibilities.

This is a series of ordinary Dirichlet FEM solves and terminal algebra, **not**
a monolithic constrained FEM formulation with floating potential unknowns.
External charge densities, nonzero dielectric-space charge, nonlinear
materials, and arbitrary impressed-flux loads are not supported here.

## Charge-constrained floating-conductor Schur solve (experimental)

`solve_floating_conductor_schur_3d` is a **discrete P1 Schur-complement
floating-conductor solver**, available through
`nova/physics/floating_conductor_schur_3d.hpp`.

Consider a linear dielectric P1 FEM system split into unconstrained interior
nodes `i` and equipotential terminal nodes `e`. Its condensed terminal
operator is

```text
C = K_ee - K_ei * inverse(K_ii) * K_ie
Q_terminal = C * V_terminal + q0
C_FF * V_F = Q_F - q0_F - C_FK * V_K
```

NOVA computes `C` implicitly through source-free unit-voltage FEM solves
and **Galerkin weak-form electrode reactions**. For impressed volumetric
free charge and/or prescribed exterior displacement `D.n`, it also solves
a loaded electrostatic problem with **every electrode at 0 V**, retaining
grounded non-electrode Dirichlet references, and evaluates the weak-form
reaction offset `q0 = (K V0 - F)_electrode`. For a source-free model it
skips that extra solve and sets `q0=0` exactly.

After solving the small affine terminal system, a final **loaded** Dirichlet
FEM solve reconstructs `V`, `E`, `D`, and field energy. The result
reports the sorted `induced_terminal_charges_coulombs` offset, terminal
voltages, **loaded** terminal-charge estimates `C V + q0`, and reconstructed
weak-form reactions. The original requested floating charges are used for
both terminal-algebra and field-reconstruction residual diagnostics.

A separate legacy direct boundary-displacement charge estimate is retained
as a mesh-sensitive comparison, not a certification of the continuum charge.
In a loaded model the direct face estimate can differ significantly from the
Galerkin `K V-F` reaction on coarse meshes. In particular the source-free
electrostatic identity `sum(V_g Q_g)=2U` must not be applied unchanged to
models with volume charge or impressed flux.

### Verification

The source-free parallel-plate regression (`epsilon=2 F/m`, unit area and
gap, `V_low=0`, `Q_high=+20 C`) returns `V_high=10 V`,
`E=(-10,0,0) V/m`, charges `[-20,+20] C`, and `U=100 J`.

An independent manufactured loaded slab has `epsilon=2 F/m`,
`rho=2 C/m^3`, `V_low=0`, `Q_high=0`. Its continuum solution is
`V(x)=x-x*x/2`, `V_high=0.5 V`, `Q_low=-2 C`, and `Q_high=0`.
The zero-terminal-voltage FEM reaction offsets are `q0=[-1,-1] C`.
The regression checks the exact potential at every structured-mesh node,
the terminal voltages, and the physically loaded charge residual.

A separate prescribed outward-displacement test applies `D.n=+1 C/m^2`
on the unit-area `y=0` side. With zero volume charge, `V_low=0`, and
`Q_high=0`, Gauss balance requires `Q_low=+1 C`. A further all-fixed
charged benchmark verifies that `estimated_charge_coulombs` includes
the affine charge offset even with no floating unknowns.

**Limits:** the solver remains linear, first-order tetrahedral, and uses
implicit Schur reduction rather than a monolithic mixed/saddle-point
formulation. The capacitance operator and its discretization are those of
the source-free P1 Galerkin model. Original non-electrode Dirichlet values
are **grounded at 0 V**, as in the capacitance basis solves, rather than
retained when nonzero; electrode groups must be disjoint, complete exterior
conductor patches. Singular or ill-conditioned floating-terminal matrices
are rejected, including disconnected gauges. This is not a fully equilibrated
charge/flux formulation. Solver and quadrature errors remain, and independent
mesh-convergence plus laboratory validation are still required.

### Analytical mesh-convergence check: graded dielectric

The `floating_schur_dielectric_capacitor_mesh_convergence` regression checks
the continuum result **independently** of the discrete Schur equations.
For a unit-area, unit-distance slab with dielectric
`epsilon(x) = (1+x) F/m`, a grounded plate at `x=0`, a floating plate
at `x=1` with `Q=+1 C`, and naturally insulated sidewalls:

```text
V_exact(x) = ln(1+x) V
V_floating = ln(2) V
C_exact = 1/ln(2) F
E_exact,x = -1/(1+x) V/m
U_exact = 0.5*ln(2) J
```

The test checks charge conservation and the FEM energy identity at multiple
mesh resolutions, then verifies that the **terminal-voltage error**, terminal
capacitance error, and potential-field error trend toward the independent
analytic solution after mesh refinement. This is a manufactured benchmark,
not experimental validation or proof of arbitrary-geometry convergence.

The earlier `reconstruct_floating_electrostatic_field_3d` based on raw
boundary flux is still available unchanged for comparisons.

## Scope and limitations

This is an experimental electrostatic Poisson interface, not a general Maxwell solver. It supports straight first-order tetrahedra, scalar or symmetric positive-definite tensor permittivity, optional per-cell material data, fixed potentials, prescribed normal displacement, and a **linear implicit Schur floating-conductor workflow**, including charge-loaded dielectric models. It does not implement an independently assembled monolithic floating-conductor constraint system, explicit surface-charge unknowns, nonlinear dielectric response, AC or transient Maxwell equations, or a locally equilibrated displacement field. Per-cell D is a constitutive element average, not an equilibrated face flux or a certified Gauss-law error estimator. This prototype is not production-certified.
