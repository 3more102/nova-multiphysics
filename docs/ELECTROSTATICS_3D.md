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

## Verification

The regression suite includes:

1. A 1 m cube parallel-plate capacitor with epsilon = 2 F/m and 10 V separation: C = 2 F and U = 100 J.
2. An anisotropic SPD dielectric tensor with an affine exact potential, exact E and D, and U = 36.5 J.
3. A one-free-DOF tetrahedron with spatial plus per-cell charge density, giving V = 1 V and U = 1/12 J.
4. The same one-free-DOF tetrahedron with cell-only charge, verifying that an omitted spatial source is treated as zero.
5. Rejection of negative permittivity, conflicting material representations and non-finite charge.

The installed SDK consumer also solves the capacitor benchmark, so the public header and exported link symbol are checked outside the source tree.

## Scope and limitations

This is an experimental electrostatic Poisson interface, not a general Maxwell solver. It supports straight first-order tetrahedra, scalar or symmetric positive-definite tensor permittivity, optional per-cell material data, fixed potentials and prescribed normal displacement. It does not yet support floating conductor potentials, surface-charge unknowns, nonlinear dielectric response, AC or transient Maxwell equations, or a locally equilibrated displacement field. Per-cell D is a constitutive element average, not an equilibrated face flux or a Gauss-law error estimator. This prototype is not production-certified.
