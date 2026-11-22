# 3D P1 diffusion residual indicator (experimental)

NOVA provides an elementwise residual indicator for the stationary scalar
problem

\[
-\nabla\cdot(k\nabla u_h)=f
\]

on conforming tetrahedral meshes. The approximation is continuous P1 and
conductivity is positive and constant within each tetrahedron (it may jump
between tetrahedra).

For each tetrahedron \(K\), the implemented indicator is

\[
\eta_K^2 =
h_K^2\int_K f^2\,dV
+\frac12\sum_{F\subset\partial K\cap\Omega}
h_F\int_F [q_h\cdot n_F]^2\,dS
+\sum_{F\subset\partial K\cap\Gamma_N^{\mathrm{selected}}}
h_F\int_F (q_h\cdot n-g_N)^2\,dS,
\qquad q_h=-k\nabla u_h.
\]

Here \(h_K\) is the longest tetrahedron edge and \(h_F\) is the longest edge
of the triangular face. The interior jump is the sum of the two outward
normal fluxes from adjacent cells. Each shared-face contribution is counted
once globally and split equally between its two incident cells.

Because a P1 gradient and the element conductivity are constant within each
cell, \(\nabla\cdot(k\nabla u_h)=0\) in the cell interior. The cell residual
is therefore just \(f\). The source-square integral is evaluated with a
five-point tensor-product Gauss rule on each Duffy-mapped reference
tetrahedron. Interior face contributions reuse the validated 3D flux-jump
diagnostic, for which the P1 normal flux jump is constant over each face.

## C++ API

Include `nova/fem/residual_error_estimator_3d.hpp` and call:

```cpp
const auto indicators = nova::fem::estimate_diffusion_residual_error_3d(
    mesh,
    nodal_values,        // one P1 value per mesh vertex
    element_conductivity,// one positive scalar k per tetrahedron
    source               // f(x,y,z) in -div(k grad(u)) = f
);
```

Pass an optional fifth argument, a vector of
`nova::fem::NeumannResidualFace3D`, to evaluate natural-boundary
flux mismatch on **explicitly selected** exterior triangles:

```cpp
std::vector<nova::fem::NeumannResidualFace3D> faces{
    {face, [](double, double, double) { return 0.0; }} // insulated face
};
const auto result = nova::fem::estimate_diffusion_residual_error_3d(
    mesh, nodal_values, element_conductivity, source, faces);
```

Outward flux is `q_h dot n = -(k grad u_h) dot n`, independent of
the supplied face winding. Each selected face adds
`h_F * ||q_h dot n - g_N||^2_L2(F)` to its owner's squared element
indicator. A positive six-point degree-four triangle rule integrates this
mismatch. Dirichlet boundaries must not be supplied; omitted natural faces
contribute no boundary residual (they are **not** assumed insulated).
Provide a zero callback explicitly for every insulated face to be included.

## Robin boundary residuals

The six-argument overload additionally accepts a vector of
`nova::fem::RobinResidualFace3D` (face, transfer coefficient `h`,
ambient field `u_infinity`). For the P1 face trace `u_h`, the implemented
Robin term is `h_F * ||q_h dot n - h(x)*(u_h(x)-u_infinity(x))||^2_L2(F)`.
The trace is interpolated from the three face vertices at each quadrature
point; it is **not** substituted by a centroid or vertex-only approximation.

`h` must be finite and nonnegative at every sampled point, and the
ambient callback must be finite. Any missing callback, duplicate face within
one set, repeated face across Neumann and Robin sets, nonexterior face, or
nonrepresentable residual is rejected. The existing four- and five-argument
overloads remain unchanged. The face report's `kind` field identifies
Neumann versus Robin residuals.

The result contains:

- `cell_residual_errors[K]`: \(h_K\|f\|_{L^2(K)}\).
- `flux_jumps`: the existing per-face and per-cell interior flux-jump report.
- `boundary_residual_errors[K]`: mismatch contribution on selected natural faces (zero in the four-argument overload).
- `boundary_residual_faces`: facewise owner, area, diameter, and L2 mismatch diagnostics.
- `element_errors[K]`: \(\eta_K\).
- `total_error`: \(\sqrt{\sum_K\eta_K^2}\).

## Scope and limitations

- P1 tetrahedral fields only; the current 3D H1 space is P1.
- Isotropic scalar conductivity, positive and constant per tetrahedron.
- Includes cell-source residuals and interior normal-flux jumps.
- Explicit Neumann (including insulated) boundary residuals are supported by the five-argument overload; the six-argument overload additionally supports Robin residuals. Natural boundaries are not automatically classified.
- Within-cell conductivity gradients, anisotropic tensor conductivity,
  temperature-dependent coefficients and nonconforming interfaces are not
  represented.
- The indicator is intended for element ranking and future marking workflows;
  it is not a guaranteed upper bound on the solution error.
- This API does not mark elements, refine a mesh or transfer a solved field.

## Verification

The automated tests check a globally affine harmonic field, the exact
constant-source indicator on a unit tetrahedron, analytical shared-face flux
jumps for a discontinuous P1 field and for discontinuous element conductivity,
Neumann matching-flux, Robin matching and mismatched flux, invalid Robin callbacks, overlapping natural face rejection, constant and spatially varying flux mismatch,
reversed face winding, duplicate/interior/missing-face rejection, and
rejection of invalid or non-finite inputs. See
[`VERIFICATION_PLAN.md`](VERIFICATION_PLAN.md).
