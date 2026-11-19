# 3D P1 tetrahedral error-norm verification

NOVA provides `compute_error_norms_3d` in
`nova/fem/error_norms_3d.hpp`. It integrates the nodal P1 FEM field
against user-supplied exact scalar and Cartesian-gradient callbacks.

The reported norms are

- `L2 = sqrt(integral (u_h-u_exact)^2 dV)`
- `H1_seminorm = sqrt(integral |grad(u_h)-grad(u_exact)|^2 dV)`
- `H1 = sqrt(L2^2 + H1_seminorm^2)`

Inputs require a validated conforming tetrahedral mesh, exactly one finite
nodal scalar value per vertex, an exact scalar function and all three
exact Cartesian derivatives. Reversed element orientation is supported.
Non-finite data and invalid sizes are rejected.

## Quadrature and limits

The integral uses the **3D Duffy map**
`r=u, s=(1-u)v, t=(1-u)(1-v)w` from the unit cube to the reference tetrahedron.
Its Jacobian is `(1-u)^2(1-v)`. A five-point Gauss-Legendre rule per
axis provides 125 integration nodes per tetrahedron. This high-order
quadrature integrates the polynomial squared error for a quadratic exact
field and its linear P1 interpolant to essentially floating-point
precision. The absolute physical determinant makes element numbering
orientation-independent.

This is a verification utility, **not** an adaptive error estimator.
It requires an analytical exact solution and may be expensive on very
large meshes. Error norms alone do not certify engineering results.

## Automated regression

`tests/test_error_norms_3d.cpp` verifies:

1. Machine-precision affine patch reproduction for both tetrahedron windings.
2. P1 interpolation of `x^2 + y^2 + z^2` over 2/4/8 box subdivisions:
   L2 error scales approximately `h^2`, H1-seminorm approximately `h`.
3. A manufactured 3D Poisson solution with `u=x^2+y^2+z^2`,
   `-div(grad u)=-6` and exact Dirichlet values on all exterior faces,
   with measured spatial refinement.
4. Invalid nodal field size, missing derivative, and non-finite data rejection.

Reproduce with `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`,
`cmake --build build --parallel`, and
`ctest --test-dir build --output-on-failure`.
