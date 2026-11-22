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
The implementation sums `-volume * grad(N_a) dot averaged_D` over all cells
and constrained electrode vertices, with the same material quadrature as the
assembled P1 stiffness. It accepts scalar, spatial or cellwise anisotropic SPD
permittivity through the regular FEM material interface.

**Required conditions:** each selected electrode group must have its complete
exterior triangular contact supplied; no vertex may belong to two electrodes;
all selected vertices must have fixed potentials and one common equipotential
per group. Nonempty free-charge callbacks, cell charge values, and prescribed
displacement-flux loads are rejected because the corresponding right-hand
side contributions are not included in this reaction diagnostic.

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

Consider the linear 3D source-free FEM stiffness matrix, split into unconstrained
interior degrees of freedom `i` and equipotential terminal degrees of freedom
`e`. Eliminating the interior unknowns yields the terminal operator

```text
C = K_ee - K_ei * inverse(K_ii) * K_ie
C_FF * V_F = Q_F - C_FK * V_K
```

NOVA obtains `C` **implicitly** through unit-terminal-voltage P1 FEM solves
and *discrete weak-form electrode reactions*, instead of constructing a dense
global nodal Schur matrix. After solving the small terminal system, a final
Dirichlet P1 solve reconstructs `V`, `E`, `D`, and stored energy.
Each returned result includes:

- The **unsymmetrized** Galerkin reaction capacitance matrix and reciprocity
  diagnostics, inferred terminal potentials, and terminal charge estimates.
- Independently integrated weak-form electrode reaction charges.
- Legacy element-averaged surface-displacement electrode estimates.
- Maximum mismatch between prescribed floating charge and reconstructed FEM
  reaction charge; maximum disagreement between weak and raw surface charges.

The parallel-plate regression (`epsilon = 2 F/m`, unit area and distance,
`V_low = 0 V`, `Q_high = +20 C`) obtains `V_high = 10 V`,
`E = (-10,0,0) V/m`, `Q = [-20,+20] C`, and `U = 100 J`.
Other regressions cover spatially varying permittivity, the discrete
`sum_g(V_g*Q_g) = 2*U` identity, an independently grounded single floating
terminal, all-fixed terminals, and singular-gauge/invalid-load rejection.

**Mathematical and engineering limits:** This is a series of standard
Galerkin FEM solves implementing an **implicit discrete Schur complement**,
not a monolithic augmented saddle-point solve with extra unknowns. A
specified independent voltage reference is required when the isolated
terminal block has a constant-voltage nullspace. It accepts a linear
permittivity representation and *no* nonempty volume-charge callback,
per-cell free charge, or prescribed-displacement faces, as those would
require an affine reaction-load term. Existing non-electrode Dirichlet
constraints are grounded at zero volts (matching the capacitance basis
solves) rather than preserving their original values. Electrode groups
must be disjoint, complete exterior conductor face patches. Solver and
quadrature tolerances limit the charge residual. The flux-versus-reaction
disagreement is **not** a validated physical error bar. Independent mesh
convergence and experimental validation remain necessary.

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

### Exact two-layer dielectric interface patch (regression)

The `floating_schur_two_layer_dielectric_exact_interface_patch` regression
checks an interface-fitted unit cube with two layers of epsilon=2 and 4 F/m
separated at x=0.5 m, insulated sidewalls, grounded x=0 electrode, and a
floating x=1 electrode with Q=+1 C. Continuity of V and normal D yields:

```text
D_x = -1 C/m^2
V(x) = x/2                     (0 <= x <= 1/2)
     = 1/4 + (x-1/2)/4        (1/2 <= x <= 1)
E_x = -1/2 V/m (left), -1/4 V/m (right)
V_right = 3/8 V
C_equivalent = 8/3 F
U = Q * V_right / 2 = 3/16 J
Q_left = -1 C; Q_right = +1 C
```

The exact potential is piecewise affine and belongs to the P1 finite-element
space on an interface-fitted mesh. The test uses 2x2x2 and 4x4x4 meshes
and repeats with both per-cell and spatial callback permittivity. It checks
individual nodal potentials and cellwise E and D; the weak-form capacitance,
floating voltage, electrode reactions, raw flux/reaction agreement, stored
energy, and the identity V*Q=2U. The test is a discretization-consistency
check, not experimental validation. It does not establish accuracy for
unfitted, curved, or nonconforming dielectric interfaces.

The earlier `reconstruct_floating_electrostatic_field_3d` based on raw
boundary flux is still available unchanged for comparisons.

## Scope and limitations

This is an experimental electrostatic Poisson interface, not a general Maxwell solver. It supports straight first-order tetrahedra, scalar or symmetric positive-definite tensor permittivity, optional per-cell material data, fixed potentials and prescribed normal displacement. It does not yet implement floating conductor boundary constraints inside the FEM solve, surface-charge unknowns, nonlinear dielectric response, AC or transient Maxwell equations, or a locally equilibrated displacement field. Per-cell D is a constitutive element average, not an equilibrated face flux or a Gauss-law error estimator. This prototype is not production-certified.
