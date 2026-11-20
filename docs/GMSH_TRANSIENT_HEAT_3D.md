# Gmsh transient 3D heat study (experimental)

NOVA's tagged 3D transient CLI combines the shared Gmsh MSH 2.2/4.1 ASCII
tetrahedral importer, the implicit P1 thermal solver, the per-cell anisotropic flux
postprocessor and the legacy VTK/PVD writers. Physical volume and exterior
surface groups can be selected by numeric ID or matching physical name. This is a reproducible
command-line study, not a GUI or an engineering certification.

## Example: two-region thermal wall

```sh
./build/nova_gmsh_transient_heat_3d examples/gmsh_diffusion_3d/two_layer_wall.msh \
  --material 11 2 --material 22 6 \
  --capacity 11 2000000 --capacity 22 3500000 \
  --initial 300 --dt 0.5 --steps 40 --scheme cn \
  --temperature 101 310 --convection 102 12 290 \
  --heat-source 22 120 --vtk-series heated_wall.pvd
```

Open `heated_wall.pvd` in ParaView. The collection contains the
Dirichlet-projected initial field at time zero and every successfully solved
time step. Each VTK frame includes nodal `temperature`, per-cell `material_id`,
effective `heat_generation` and element-averaged `heat_flux`. The exporter
uses one FEM assembly and one transient solve sequence; it does not replay or
re-solve the simulation to produce frames. The CLI disables retention of
intermediate nodal solution vectors in every output mode; PVD frames are
written synchronously by the solver observer.

For a single final result, use `--vtk heated_wall_final.vtk` instead of
`--vtk-series heated_wall.pvd`.

## Inspect mesh quality before running heat studies

The transient Gmsh CLI uses the same validated P1 tetrahedral mesh-quality
metric as the steady 3D CLI. Inspect a mesh without specifying materials or
running a thermal solve:

```sh
./build/nova_gmsh_transient_heat_3d examples/gmsh_diffusion_3d/two_layer_wall.msh --check-mesh
```

The report includes total volume, minimum mean ratio (0 to 1), worst-element
index and minimum dihedral angle. Add `--quality-report` to a normal thermal
study to print these diagnostics. Add `--min-mean-ratio q` with
`0 < q <= 1` to print the report and reject the mesh **before FEM assembly or
any VTK/PVD result writes** if its worst tetrahedron has mean ratio below
`q`. The quality-only `--check-mesh` mode accepts the same threshold but
cannot be combined with study settings.

This is a geometric screening threshold, not a guarantee of convergence,
accuracy or engineering suitability. Choose a tolerance appropriate to the
meshing workflow and verify the solved fields independently.

## Inputs and conventions

- `--material VOL k`: finite positive isotropic conductivity (W/m/K) for a
  tetrahedral physical volume tag.
- `--material-tensor VOL xx yy zz xy xz yz`: symmetric positive-definite
  conductivity tensor in global Cartesian coordinates. Isotropic and tensor
  materials can be mixed across volumes.
- `--capacity VOL C`: positive volumetric heat capacity
  `C = rho * cp`, in J/(m³ K). Every tetrahedral physical volume tag needs
  exactly one material and one capacity.
- `--initial T0`: uniform initial nodal temperature (K). Dirichlet face
  temperatures are projected onto constrained nodes at time zero.
- `--source Q`: uniform volumetric heating (W/m³), default zero.
  `--heat-source VOL Q` adds a signed source density to the selected volume.
- `--temperature SURF T`: fixed temperature (K) on all exterior triangles
  carrying the selected physical surface tag.
- `--flux SURF q_out`: prescribed outward heat flux (W/m²); positive values
  remove heat from the domain.
- `--convection SURF h T_inf`: Robin boundary with outward flux
  `q_out = h (T - T_inf)`, where `h >= 0` is in W/(m² K).
- `--dt seconds`, `--steps N`: constant positive time step and step count.
- `--scheme be|cn`: Backward Euler (default) or Crank–Nicolson.
- `--vtk output.vtk`: optional final-state VTK file. It is written through a
  temporary sibling directory and renamed after a complete file write.
- `--vtk-series output.pvd`: write the projected initial state and all
  completed steps as timestamped VTK frames plus a PVD collection.

The importer accepts first-order four-node tetrahedra and triangular exterior
facets from Gmsh **MSH 2.2 or 4.1 ASCII**. MSH 4.1 physical groups are mapped
from geometric entities, including named groups, sparse identifiers and
multiple node/element blocks. Unassigned exterior faces are insulated.
Duplicate or conflicting boundary assignments, missing material/capacity
coverage, unknown physical tags, invalid tensor conductivity, nonpositive
capacity, non-finite values and output-path collisions are rejected.
`--vtk` and `--vtk-series` are mutually exclusive. The PVD destination must
have a safe ASCII basename and an existing parent directory. Existing PVD or
frame files are not overwritten. Individual frames and the collection are
staged and renamed after successful writes; if a later solve fails, already
published frames remain without a PVD collection.

## Diagnostics and verification

The CLI reports the final temperature range and the maximum absolute
per-step discrete FEM energy-budget residual. The residual includes
Dirichlet reactions and is an algebraic conservation diagnostic; it is not
an independent measurement of physical boundary heat flux.

The `nova_gmsh_transient_heat_3d_verify` CTest performs an end-to-end
two-material insulated-wall study with spatially uniform capacity `C=2`
and source `Q=4`, for which `T(t) = 300 + 2t` is exact. It parses the PVD XML
and every generated VTK frame, checks all physical times and nodal
temperatures, and verifies the cell output fields. It also exercises a mixed
scalar/tensor material model with Dirichlet and convection boundaries,
final-only VTK output, invalid groups and material inputs, and overwrite
protection.

## Scope limits

This is a fixed-step, time-independent-coefficient P1 prototype. It does not
yet provide adaptive 3D time stepping, time-varying material/boundary data,
temperature-dependent conductivity, phase change, radiation, moving meshes,
restart checkpoints, higher-order tetrahedral elements or certified
engineering validation. MSH 4.1 binary files, partitioned/periodic metadata,
unsupported sections, high-order elements and entities with multiple physical
groups are rejected. The current MSH 4.1 adapter buffers and normalizes the
input before topology validation, so its peak memory is higher than a future
streaming native reader.
