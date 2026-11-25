# 2D scalar waveguide eigenmodes (experimental)

## Scope

NOVA provides a P1 triangular finite-element eigenmode reference for a slab waveguide propagating along the global x-axis, with a transverse y-z cross-section. The generic coordinates of `TriangleMesh2D` parameterize that cross-section. Materials are isotropic, lossless, and piecewise-constant, with positive real relative permittivity and permeability. It solves the two scalar reductions:

- `TE_Ez`: the out-of-plane electric field component (E_z).
- `TM_Hz`: the out-of-plane magnetic field component (H_z).

For this propagation-axis convention, TE_Ez has a transverse z-directed electric field and TM_Hz has a transverse z-directed magnetic field. These labels should not be confused with unrelated plane-of-incidence conventions. This is **not** a full-vector hybrid-mode solver.

## Discrete formulation

Let (k_0=2\pi/\lambda_0), (\beta=k_0 n_\mathrm{eff}), and define

`K_a[i,j] = integral(a * grad(phi_i) dot grad(phi_j))`
and
`M_a[i,j] = integral(a * phi_i * phi_j)`.

The generalized symmetric eigenproblems are

- TE_Ez: `(k0^2 M_epsilon - K_(1/mu)) e = beta^2 M_(1/mu) e`.
- TM_Hz: `(k0^2 M_mu - K_(1/epsilon)) h = beta^2 M_(1/epsilon) h`.

The implementation uses exact P1 triangle stiffness and consistent mass integrals for cellwise-constant material coefficients. The sparse generalized Krylov/Rayleigh-Ritz solver applies the generalized operator through Jacobi-preconditioned conjugate-gradient mass solves and two-pass M-inner-product orthogonalization. It forms a bounded dense projection, reconstructs Ritz vectors in the original sparse space, and checks their true generalized residuals. Eigenvectors are normalized by the corresponding weighted mass matrix. Returned candidates are ordered by descending effective index and carry a weighted core-confinement factor (modal mass in cells whose refractive index exceeds the declared cladding index).

## Boundary treatment

The computational domain is finite. If `dirichlet_boundary_nodes` is empty, all exterior vertices receive homogeneous scalar Dirichlet values. If the list is supplied, every listed node must be an exterior vertex; unlisted exterior edges receive the natural homogeneous Neumann condition. The latter supports symmetry planes and invariant directions in slab benchmarks. It is not an exact open-radiation boundary condition.

Only eigenpairs satisfying `cladding_index < n_eff < max(n_material)` are returned. This light-line filter alone cannot distinguish every physical guided mode from a finite-window discretization of the cladding continuum. Inspect `core_confinement_factor` and repeat the calculation with a larger cladding window before interpreting near-cutoff candidates as physical guided modes.

## Numerical limits and exclusions

- P1 triangular H1 scalar fields only.
- Real, positive, isotropic, cellwise-constant (\epsilon_r) and (\mu_r).
- Sparse generalized Krylov/Rayleigh-Ritz solve with a dense projected subspace limited to 512 dimensions; the full sparse system is not subject to the old 512-DOF dense limit.
- Finite computational window with homogeneous Dirichlet and natural homogeneous Neumann boundaries.
- No complex/lossy materials, anisotropy, PML, ports, radiation boundary operators, curved geometry, higher-order elements, or full-vector hybrid modes.
- The Krylov method is currently unrestarted. Its subspace is bounded by the configured projection dimension (default 64); insufficient subspace convergence is reported through `converged=false`. Jacobi-preconditioned CG mass solves may require more iterations for strongly graded meshes. The confinement factor is a diagnostic, not a certified modal classification or an error bound.

This is an experimental numerical reference and is not qualified for engineering design or commercial release.

## Verification

`tests/test_waveguide_modes_2d.cpp` compares the fundamental TE_Ez and TM_Hz modes of a symmetric slab against an independently bisected analytical dispersion relation. The test checks mesh refinement against an independent analytical slab reference, material-weighted TE/TM operators, mode ordering, true residuals, boundary constraints, natural side boundaries, triangle-orientation invariance, homogeneous-cladding no-mode behavior, and invalid inputs.

Build and run:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The `nova_waveguide_modes_2d` example runs both scalar polarizations on the same slab cross-section.
