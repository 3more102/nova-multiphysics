# Spatially varying impressed Maxwell sources

The installed header nova/physics/maxwell_spatial_source_3d.hpp extends
the driven 3D H(curl) solver without modifying its assembly or solve APIs.
The callback receives the cell index and physical SI coordinates and returns
F in curl(mu_r^-1 curl E)-k0^2 epsilon_r E=F. For impressed electric current
density J with exp(-i omega t), F=i omega mu0 J.

Example:
```cpp
const auto load = nova::physics::assemble_maxwell_spatial_load_3d(
    mesh, [](std::size_t, const std::array<double,3>& p) {
        return nova::physics::ComplexVector3D{std::sin(p[2]),0.0,0.0};
    }, system);
const auto solution = nova::physics::solve_maxwell_3d(system,load);
```

The Duffy map (r,s,t)=(u,(1-u)v,(1-u)(1-v)h) has Jacobian
(1-u)^2(1-v). Four-point Gauss-Legendre on each coordinate of the unit
cube is exact up to degree seven per coordinate. Consequently the
tetrahedral rule integrates total polynomial degree five exactly in
exact arithmetic, including a quartic vector source times the affine
Nedelec basis. Physical Jacobian is 6 times tetrahedron volume.
64 callback evaluations per cell; the existing centroid API is
preferable for constant sources. Returns full edge loads, and the
existing solver handles PEC reduction.

Tests: heterogeneous piecewise constant equivalence, unit-tetrahedron
quartic monomial exact integrals and orientation invariance, missing
callback, nonfinite source and mismatched edge topology. Polynomial
identity is a verification target, not a measured full-wave convergence
result. No adaptive quadrature, singular source, point dipole, surface
current, port, or PML implementation. For discontinuities, align the
mesh to the interface or use a suitable integration method. Smooth
nonpolynomial sources may require mesh refinement to resolve variation.
