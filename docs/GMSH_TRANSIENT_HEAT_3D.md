# Imported Gmsh transient 3D heat study (experimental)

NOVA's tagged transient CLI combines the validated Gmsh MSH 2.2 ASCII
first-order tetrahedral importer, the P1 transient heat solver and the
ParaView PVD/legacy VTK writers. It is a reproducible command-line study,
not a CAD front end or a certified engineering workflow.

## Example: two-region thermal wall

~~~sh
./build/nova_gmsh_transient_heat_3d examples/gmsh_diffusion_3d/two_layer_wall.msh \
  --material 11 2 --material 22 6 \
  --capacity 11 2000000 --capacity 22 3500000 \
  --initial 300 --start-time 0 --dt 0.5 --steps 40 --scheme cn \
  --temperature 101 310 --convection 102 12 290 \
  --vtk-series heated_wall.pvd
~~~

Open heated_wall.pvd in ParaView. The collection contains frame zero after
initial Dirichlet projection and one frame after each successfully converged
time step. The solver assembles the fixed operators once; the observer writes
the actual solved state and does not replay or re-solve the trajectory.
retain_step_solutions=false prevents retaining an O(steps × nodes) array of
intermediate nodal temperatures.

## Input options and conventions

- --material VOL k: positive isotropic conductivity (W/(m K)) for a
  tetrahedral physical volume tag.
- --material-tensor VOL xx yy zz xy xz yz: symmetric positive-definite
  conductivity tensor in global Cartesian coordinates. Scalar and tensor
  regions may be mixed.
- --capacity VOL C: positive volumetric heat capacity (J/(m³ K)); every
  tetrahedral volume tag must have exactly one conductivity and capacity.
- --heat-source VOL Q: additive, signed region source in W/m³.
  --source Q adds a uniform source to all cells.
- --initial T0: uniform initial nodal temperature (K).
- --start-time seconds: physical time associated with the initial field
  (default 0).
- --dt seconds, --steps N: positive fixed time step and step count.
- --scheme be|cn: Backward Euler or Crank–Nicolson.
- --temperature SURF T: fixed temperature (K) on an exterior physical
  surface group.
- --flux SURF q_out: prescribed outward heat flux (W/m²); positive is
  heat leaving the solid.
- --convection SURF h T_inf: outward Robin law
  q_out = h(T - T_inf), with nonnegative h in W/(m² K).
- --vtk output.vtk: optional final-state VTK export.
- --vtk-series output.pvd: initial plus every completed state as a ParaView
  time series. --vtk and --vtk-series are mutually exclusive.

VOL and SURF accept numeric physical tags or exact names from Gmsh
$PhysicalNames, with dimension-specific lookup. Unassigned exterior faces
are insulated. Duplicate/conflicting surface assignments, unknown tags,
missing region material/capacity, nonpositive capacities, non-SPD tensors and
nonfinite numeric inputs are rejected.

Each VTK frame contains nodal temperature, cell material_id,
heat_generation and element-averaged physical heat_flux. PVD times are
the physical solver times, not frame indices. The collection manifest is
published only after every frame has been written; deterministic output
names are checked for collisions before solving. Existing output files are
not overwritten.

## Verification

nova_gmsh_transient_heat_3d_pvd_verify independently parses the generated
PVD XML and every VTK frame using only the Python standard library. It checks:

- A two-material insulated heating case with C=2, Q=4, and
  T(t)=300+2(t-t0), including the exact frame times, cell tags,
  heat-generation field and zero heat flux.
- Equality of the final PVD frame and a separately exported final VTK field.
- A mixed scalar/tensor material case with Dirichlet and convection groups,
  regional heat generation and a closed discrete energy budget.
- A named physical-group case using the unit-tetrahedron fixture.
- Rejection of malformed output paths, duplicate output modes, missing
  materials, invalid capacities and unknown boundary tags.

The reported energy residual is the discrete FEM heat-budget closure and
includes the Dirichlet reaction. It is not an independent measurement of
physical boundary heat flux.

## Scope limits

This path currently accepts Gmsh MSH 2.2 ASCII, first-order tetrahedra and
time-independent material/source/boundary data with a fixed time step. It
does not provide transient MSH 4.1 import, adaptive time stepping, nonlinear
or temperature-dependent materials, radiation, moving meshes, restart CLI,
CAD geometry editing, experimental validation or engineering certification.
