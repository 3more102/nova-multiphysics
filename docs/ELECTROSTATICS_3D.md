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

## Scope and limitations

This is an experimental electrostatic Poisson interface, not a general Maxwell solver. It supports straight first-order tetrahedra, scalar or symmetric positive-definite tensor permittivity, optional per-cell material data, fixed potentials and prescribed normal displacement. It does not yet support floating conductor potentials, surface-charge unknowns, nonlinear dielectric response, AC or transient Maxwell equations, or a locally equilibrated displacement field. Per-cell D is a constitutive element average, not an equilibrated face flux or a Gauss-law error estimator. This prototype is not production-certified.
