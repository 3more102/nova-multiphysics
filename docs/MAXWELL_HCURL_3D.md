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

No eigenvalue solver, eigenmode/spurious-mode validation, dispersive material database, nonisotropic tensors, PML, port, incident wave boundary, selective PEC, adaptive H(curl), or nonlinear optical coupling. Lossless resonant systems may be singular; a converged algebraic solve is not evidence of physical convergence. The unpreconditioned/restarted GMRES approach is not a scalable Maxwell solver; no performance claims. No cavity or waveguide accuracy claim is made. Meshes with every edge constrained are rejected. Piecewise constant isotropic passive epsilon (including negative real parts) and positive real mu are accepted; near-zero or strongly indefinite systems need additional verification.

Next gates: non-polynomial manufactured field and h-convergence, cavity eigenfrequency with gradient-space filtering, selective boundary conditions, and scattering/PML verification against existing optical reference models.
