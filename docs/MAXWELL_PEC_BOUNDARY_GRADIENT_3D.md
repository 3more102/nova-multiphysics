# PEC cavity boundary-component gradient projection

For lowest-order Nédélec edges on a mesh with full PEC conditions, the edge
circulations of a piecewise linear scalar potential are endpoint differences.
A PEC-compatible scalar potential is constant on each connected boundary
surface. Distinct surfaces can have different constants. Fixing all of them
to zero omits static gradients, e.g. the field between an enclosed conducting
surface and the outer cavity wall.

The sparse solver now assigns one independent scalar potential to each
interior node and to each PEC surface except one reference surface per
connected mesh component. It constructs the sparse incidence G and uses
the existing mass-orthogonal projection
P = I - G (G^T M G)^-1 G^T M.
The reference fixes the constant gauge without imposing an extra physical
constraint. `boundary_component_gradient_dimension` reports the number
of added surface-potential columns; `gradient_nullspace_dimension`
includes both interior and boundary-component columns.

The hollow-cube regression independently forms the exact discrete static
field from unit inner-surface and zero outer-surface potentials, checks
K g = 0 and g^T M g > 0, compares positive sparse modes with the existing
dense reference, and checks mass orthogonality to g.

This is a discrete scalar-gradient construction, not a general cohomology
or topology certification. Mesh connectivity is vertex/edge based; meshes
with touching nonmanifold surfaces require separate qualification. Additional
unprojected zero modes are still reported by the existing Ritz filter.
The new regression is pending hosted CI at submission.
