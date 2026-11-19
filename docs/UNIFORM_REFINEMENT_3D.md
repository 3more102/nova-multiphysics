# Conforming uniform tetrahedral refinement (P1)

NOVA refines a validated, conforming first-order tetrahedral mesh into eight
children per parent without introducing hanging nodes. This is a global
uniform refinement tool, **not** error-driven adaptive mesh refinement.

## Mathematical construction

1. Create a globally shared midpoint for every unique original mesh edge.
   Midpoint numbering uses lexicographically sorted coarse vertex IDs, so
   the numbering does not change with tetrahedron insertion order.
2. Split every parent tetrahedron into four corner tetrahedra plus one
   central octahedron.
3. Split each octahedron into four tetrahedra along its shortest interior
   diagonal among three choices. The chosen diagonal affects internal mesh
   quality but never the triangulation of a parent exterior face.
4. For every supplied outward exterior triangle, generate four outward
   triangles and report the original parent boundary-face index.
5. Validate the complete child mesh using the existing 3D manifold checks.

For a mesh with `V` vertices, `E` unique edges and `T` tetrahedra,
the refined mesh has `V+E` vertices and `8T` tetrahedra. A supplied exterior
face set of `F` triangles becomes `4F` triangles. The implementation
rejects invalid, duplicated, reversed, or non-exterior supplied boundary
faces. Unspecified exterior faces are still present geometrically but are
not included in the returned `boundary_faces` metadata.

Each child stores its parent cell ID. An input field that is constant
per parent cell is copied to every child. A coarse continuous P1
field is interpolated to the child mesh at original vertices and edge
midpoints, which exactly reproduces global affine fields.

## Usage

```cpp
#include "nova/fem/mesh_refinement_3d.hpp"

const auto refined = nova::fem::refine_tetrahedra_uniformly_3d(
    coarse_mesh, outward_boundary_faces
);
const auto fine_temperature = nova::fem::prolongate_p1_uniform_3d(
    refined, coarse_temperature
);
const auto fine_conductivity =
    nova::fem::prolongate_cellwise_constant_uniform_3d(
        refined, coarse_cell_conductivity
    );
// Rebuild a physics problem using refined.mesh and boundary faces; do not
// reuse the original mesh's connectivity or node indices.
```

The `parent_boundary_faces` array maps each child face to the supplied
boundary array, allowing caller-managed physical tags and boundary data to
be replicated without silently changing their meaning. The
`parent_elements` array maps child cells to original cells.

## Verification and limitations

C++ regression tests cover a single parent tetrahedron, a pair of
face-sharing tetrahedra, two successive refinement levels, exact
parent volume conservation, exterior-face set equality and orientation,
affine P1 interpolation, element-constant transfer and invalid inputs.
These are implementation checks, not industrial-grade solution accuracy
or mesh-quality guarantees.

Uniform refinement multiplies the tetrahedron count by eight per level;
memory and compute costs therefore grow rapidly. The algorithm does not
generate local/graded adaptive meshes, transfer higher-order DOFs,
automatically remap physical condition objects, or optimize distorted
elements. A subsequent step should use this API for mesh-convergence
benchmarks with independent manufactured solutions.
