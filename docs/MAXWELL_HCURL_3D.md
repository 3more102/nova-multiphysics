# Driven 3D Maxwell FEM foundation

This experimental C++20 module uses the existing tetrahedral mesh, physical P1 gradients, real CSR and restarted GMRES. No new third-party dependency is introduced. It does not alter the existing optics reference branches.

## Formulation

SI coordinates and vacuum wavenumber k0 = omega/c; phasor convention exp(-i omega t). On a piecewise homogeneous isotropic material solve

curl(mu_r^-1 curl E) - k0^2 epsilon_r E = F.

The weak form integrates mu_r^-1 curl(w) dot curl(E) - k0^2 epsilon_r w dot E against real H(curl) basis functions. The natural boundary is homogeneous n cross (mu_r^-1 curl E) = 0. Optional homogeneous PEC on all exterior faces eliminates their edge unknowns. For an impressed current density J, F = i omega mu0 J under the stated convention.

The first-family lowest-order basis on each affine tetrahedron is w_ij = s_ij (lambda_i grad(lambda_j) - lambda_j grad(lambda_i)), with curl(w_ij) = 2 s_ij grad(lambda_i) cross grad(lambda_j). The sign makes every global edge point from the smaller node index to the larger. Its line integral on that edge is one and on all other edges is zero. The curl-curl integral is V times the curl dot product, and the mass integral uses integral(lambda_i lambda_j)=V/10 for i=j, V/20 otherwise. Constant vector sources use exact centroid integration. Complex A is solved through the equivalent real CSR block [Re A, -Im A; Im A, Re A].

## Verification gates

- Exact unit-tetrahedron mass and curl-curl entries, discrete gradient nullspace and oriented edge circulation.
- Global assembly invariant under a local vertex permutation.
- Manufactured constant complex electric field on 6 and 48 tetrahedra with F=-k0^2 epsilon_r E, natural boundary, and nonzero material loss. Exact edge circulation is E dot (x_b-x_a); the test requires maximum absolute error below 1e-9 with GMRES relative tolerance 1e-12. This is a polynomial reproduction test, not a measured asymptotic convergence order.
- Exterior PEC edge elimination on a six-tetrahedron cube and zero-load recovery.
- Rejection of active media and invalid frequency.

Run the full existing suite with:
```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNOVA_BUILD_TESTS=ON
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

For manufactured errors, run ctest with -V -R nova_unit_tests.

## Limits

No eigenvalue solver, eigenmode/spurious-mode validation, dispersive material database, nonisotropic tensors, PML, port, incident wave boundary, adaptive H(curl), or nonlinear optical coupling. Lossless resonant systems may be singular; a converged algebraic solve is not evidence of physical convergence. The unpreconditioned/restarted GMRES approach is not a scalable Maxwell solver; no performance claims. No cavity or waveguide accuracy claim is made. Meshes with every edge constrained are rejected. Piecewise constant isotropic passive epsilon (including negative real parts) and positive real mu are accepted; near-zero or strongly indefinite systems need additional verification.

Next gates: non-polynomial manufactured field and h-convergence, cavity eigenfrequency with gradient-space filtering, and scattering/PML verification against existing optical reference models.

## Field reconstruction and optical heating

For peak phasor E, reconstruct H = curl(E)/(i omega mu0 mu_r). The time-averaged power flux is Re(E cross conjugate(H))/2 in W/m^2. For nonmagnetic-loss passive isotropic materials, the dielectric absorption density is q = omega epsilon0 Im(epsilon_r) |E|^2 / 2 in W/m^3. The API returns per-cell absorbed power in watts using a positive four-point degree-two rule, exact for the affine edge field. It can serve as a conservative source for a later thermal coupling workflow, but no thermal feedback is implemented here.

Added analytical verification: E=(-y,x,i) on the unit tetrahedron; curl, magnetic field, loss density and integrated loss are checked independently. The exact loss integral is omega epsilon0 Im(epsilon_r) * (1+2/10)/12. A second complex affine field checks nonzero Poynting flux and lossless behavior. These are field reconstruction and integration tests, not an optical scattering energy-balance or coupled thermal benchmark.

## Smooth PEC manufactured h-convergence

On the unit cube take E=(0,sin(pi*x) sin(pi*z),0), with homogeneous exterior PEC. This field is divergence-free, has curl-curl E=2*pi^2 E, and satisfies the tangential boundary condition. With k0=2.3 and epsilon_r=2+0.6i, F=(2*pi^2-k0^2 epsilon_r)E. The implemented cellwise-constant source samples F at each tetrahedron centroid. A five-point tensor Gauss rule on a Duffy-transformed tetrahedron estimates the electric-field L2 error.

The independent NumPy/SciPy reference was executed using direct sparse complex solution, high-order numerical mass integration (instead of the closed form), and the same centroid-source approximation. Results:
| Tetrahedra | Free edges | Relative L2 error | Relative algebraic residual |
| ---: | ---: | ---: | ---: |
| 48 | 26 | 0.618245911 | 1.64e-15 |
| 384 | 316 | 0.318023797 | 7.63e-15 |
| 3072 | 3032 | 0.159994250 | 5.17e-14 |

Observed orders are 0.959 and 0.991. These are independent reference results, NOT measured C++ results. The C++ regression test requires the first two errors to agree to 1e-6 absolute and the observed order to exceed 0.8, with GMRES relative residual tolerance 1e-11. The high-order quadrature is a numerical reference for the error, not an exact integral of the sine field. The finest grid still has approximately 16% error, not a production optical accuracy claim.

The optional reference script requires NumPy and SciPy only for verification, not for the C++ build or SDK:
```sh
python verification/verify_maxwell_hcurl_reference.py
```

## Selective homogeneous PEC

`assemble_maxwell_3d_selective_pec(mesh, materials, k0, pec_faces)`
accepts a list of exterior triangles identified by their three mesh node
indices. Every edge belonging to a selected face is constrained to zero
circulation, including edges shared with unselected exterior faces.
Unselected faces retain the homogeneous natural curl condition, subject
to those shared-edge constraints. Face orientation and list order do
not matter. Duplicate triangles, internal interfaces, degenerate faces
and nonexistent faces are rejected. Empty and full selections reproduce
the existing natural and all-exterior PEC assembly respectively.

Regression gates include a complex constant E normal to a selected face
of a unit tetrahedron, with F=-k0^2 epsilon_r E; exact edge circulation
and zero constrained values are checked. Additional tests cover legacy
matrix equivalence, reversed orientations and invalid face selection.
These are proposed regression gates, not a new measured convergence study.
No nonhomogeneous tangential excitation, impedance condition, ports or
absorbing boundary is provided by this API.


## Exterior Poynting power diagnostic

`integrate_maxwell_boundary_flux_3d(mesh, materials, k0, system, edge_field)`
returns one signed time-averaged power in watts per exterior triangle,
with sorted node IDs and the owning tetrahedron index. Positive values
leave the domain. Sum selected face values to measure net power through
a mesh-aligned exterior surface. The area vector is oriented away from
the owner tetrahedron; local vertex ordering does not affect the result.

For lowest-order affine tetrahedral H(curl) fields, E is affine and
H is constant in each cell. Therefore S = Re(E cross conjugate(H))/2
is affine and centroid integration of S dot n is exact on a planar face.
The method excludes internal faces and rejects nonmanifold face incidence,
incompatible edge ordering, nonfinite coefficients and invalid materials.
The regression fixtures check signed analytic flux, permutation invariance,
zero net flux of a constant S, shared-face exclusion and invalid input.

This is a one-sided reconstructed field diagnostic, not an equilibrated
flux. Do not assume an exact global power balance for arbitrary discrete
fields. No absorbing boundary, incident/scattered separation, waveguide
port normalization, S-parameter extraction or PML is implemented by this
API. Full C++ build and test execution are required before merging.
