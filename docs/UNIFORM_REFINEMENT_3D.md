# Conforming uniform tetrahedral refinement (experimental)

`nova/fem/mesh_refinement_3d.hpp` defines `refine_tetrahedra_uniform_3d`
and `prolongate_tetrahedral_p1_3d`.

Every conforming straight-sided P1 tetrahedron is split into **eight**
children using the six midpoints of its original edges: four corner tetrahedra
and four tetrahedra obtained by splitting the interior octahedron along an
opposite-edge midpoint diagonal. An edge-key map ensures adjacent original
cells reuse a single midpoint, making all shared face subdivisions conforming.

`UniformTetrahedralRefinement3D` contains:
- `mesh`: the validated refined tetrahedral mesh.
- `original_node_count` and `midpoint_parent_edges`: a deterministic map
  from coarse vertices and edges to fine P1 nodal degrees of freedom.
- `parent_elements`: one original tetrahedron index per refined element,
  usable to transfer elementwise material IDs or coefficients.
- `boundary_faces` and `parent_boundary_faces`: four outward-oriented child
  triangles for each coarse exterior triangle, indexed in the order returned
  by `extract_boundary_faces_3d(coarse)`. This preserves parent surface
  attribution for future physical boundary-tag transfer.

`prolongate_tetrahedral_p1_3d` exactly interpolates the coarse finite-element
field at each new edge midpoint, sharing values across element interfaces.
In particular, arbitrary affine scalar fields are reproduced exactly at
every refined vertex (up to roundoff).

Example:

```cpp
#include "nova/fem/mesh_refinement_3d.hpp"
auto coarse = nova::fem::make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
auto refined = nova::fem::refine_tetrahedra_uniform_3d(coarse);
std::vector<double> temperatures(coarse.nodes.size(), 300.0);
auto refined_temperatures =
    nova::fem::prolongate_tetrahedral_p1_3d(refined, temperatures);
```

The regression suite verifies eight equal-volume children on both tetrahedron
orientations, face conformity across two tetrahedra, total-volume conservation,
coarse-to-fine provenance and exact affine P1 interpolation.

**Scope:** this is uniform 1-to-8 tetrahedral refinement, not error-indicator
driven adaptivity. It does not remesh curved CAD surfaces, introduce higher
order elements, automatically rebind physics boundary conditions, or transfer
solver matrices. Transferred physical tags are represented by parent index
mappings rather than embedded in the bare mesh.
