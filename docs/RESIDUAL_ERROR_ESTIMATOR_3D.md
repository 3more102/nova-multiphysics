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
h_F\int_F [q_h\cdot n_F]^2\,dS,
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

The result contains:

- `cell_residual_errors[K]`: \(h_K\|f\|_{L^2(K)}\).
- `flux_jumps`: the existing per-face and per-cell interior flux-jump report.
- `element_errors[K]`: \(\eta_K\).
- `total_error`: \(\sqrt{\sum_K\eta_K^2}\).

## Scope and limitations

- P1 tetrahedral fields only; the current 3D H1 space is P1.
- Isotropic scalar conductivity, positive and constant per tetrahedron.
- Includes cell-source residuals and interior normal-flux jumps.
- Natural-boundary residuals (Neumann or Robin mismatch) are omitted.
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
and rejection of invalid or non-finite inputs. See
[`VERIFICATION_PLAN.md`](VERIFICATION_PLAN.md).
