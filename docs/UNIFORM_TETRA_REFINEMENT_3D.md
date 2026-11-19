# Uniform tetrahedral refinement (experimental)

NOVA provides conforming global 1-to-8 refinement of straight-sided P1
tetrahedral meshes. This operation is independent of the steady and transient
3D heat solvers; it produces new mesh geometry and transfer maps, but does not
silently alter the PDE coefficients or boundary conditions.

## C++20 usage

```cpp
#include "nova/fem/mesh_refinement_3d.hpp"

auto box = nova::fem::make_box_mesh_3d(
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 1, 1
);
auto refined = nova::fem::refine_tetrahedra_uniformly_3d(
    box.mesh, box.boundary_faces
);
std::vector<double> temperatures(box.mesh.nodes.size(), 300.0);
auto refined_temperatures =
    nova::fem::prolongate_p1_uniform_3d(refined, temperatures);
```

Each unique coarse edge gets **one shared midpoint**, so neighboring
tetrahedra do not create hanging nodes. Four corner tetrahedra and four
tetrahedra filling the central octahedron replace each parent. Among the three
possible central octahedron diagonals, the shortest geometric diagonal is
chosen; this does **not** prove optimal aspect ratio or minimum angle.

- `parent_elements` identifies the coarse tetrahedron for every fine cell.
- `midpoint_parent_edges` allows exact interpolation of P1 nodal fields
  for affine functions.
- `parent_boundary_faces` identifies which supplied outward-oriented coarse
  exterior face generated each of four fine boundary triangles.
- `prolongate_cellwise_constant_uniform_3d` repeats finite per-element
  scalar values, e.g. conductivity, volumetric heat capacity or source density,
  across each parent's children.

The source mesh must pass `TetrahedralMesh3D::validate()`. Supplied boundary
faces must be distinct, outward-oriented, and exterior. A partial exterior
face list is allowed. Refined connectivity is validated before return.

## Verification expectations

The unit tests check multi-cell conforming topology, preserved domain
volume, exact oriented exterior face subdivision, parent maps, affine P1
interpolation across two refinement generations, preservation of integrated
piecewise-constant source strength, reversed parent orientation, and invalid
boundary/transfer inputs.

Run the full project suite:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Explicit limitations

- Uniform global refinement only: no local adaptivity, hanging-node
  constraints, coarsening, or dynamic repartitioning.
- P1 geometry/nodal transfer and scalar cellwise transfer only; no P2
  transfer, finite-element solution error estimate, automatic material
  database transfer, or boundary-condition callback remapping.
- Midpoints assume straight edges; curved boundaries require a separate
  geometry projection step.
- This is numerical infrastructure under development, **not**
  independently certified engineering meshing software.
