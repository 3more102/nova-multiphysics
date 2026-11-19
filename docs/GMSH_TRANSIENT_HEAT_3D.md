# Gmsh transient 3D heat study (experimental)

NOVA's tagged 3D transient CLI combines the verified Gmsh MSH 2.2 tetrahedral
importer, the implicit P1 thermal solver and the legacy VTK conduction writer.
It accepts physical volume and exterior surface **IDs** or names. It is a
reproducible command-line study, not a GUI, a PVD time series or engineering
certification.

## Example: two-region thermal wall

```sh
./build/nova_gmsh_transient_heat_3d examples/gmsh_diffusion_3d/two_layer_wall.msh \
  --material 11 2 --material 22 6 \
  --capacity 11 2000000 --capacity 22 3500000 \
  --initial 300 --dt 0.5 --steps 40 --scheme cn \
  --temperature 101 310 --convection 102 12 290 \
  --heat-source 22 120 --vtk heated_wall_final.vtk
```

The output includes nodal `temperature`, per-cell `material_id`, effective
`heat_generation` and `heat_flux` fields readable by ParaView. The command
prints the minimum and maximum final temperatures and the maximum absolute
*discrete FEM energy-budget residual* over all steps. The energy residual
includes Dirichlet reactions; it is **not** an independent measurement of
physical boundary heat flux.

## Inputs and conventions

- `--material VOL k`: one finite positive isotropic conductivity (W/m/K)
  for every tetrahedral physical volume tag.
- `--material-tensor VOL xx yy zz xy xz yz`: symmetric, positive-definite
  conductivity in global Cartesian coordinates. Isotropic volumes can be
  mixed with anisotropic volumes.
- `--capacity VOL C`: finite positive volumetric heat capacity
  `C = rho * cp`, in J/(m³ K), required for **every** volume tag.
- `--initial T0`: uniform nodal initial temperature in K. Prescribed
  Dirichlet temperature takes priority on constrained nodes at `t=0`.
- `--source Q`: uniform volumetric heating in W/m³ (default 0).
  `--heat-source VOL Q` adds a signed source on the specified volume.
- `--temperature SURF T`: fixed exterior temperature (K).
- `--flux SURF q_out`: positive heat flow **outward** (W/m²).
- `--convection SURF h T_inf`: outward Robin heat flux
  `q_out = h (T − T_inf)`, with `h >= 0` in W/(m² K).
- `--dt seconds`, `--steps N`: constant positive time step and step count.
- `--scheme be|cn`: Backward Euler (default) or Crank–Nicolson.
- `--vtk output.vtk`: optional **final** solution, written via a temporary
  sibling directory and renamed after a complete file write.

Each assigned physical surface tag applies to all corresponding exterior
triangles. Missing boundary tags are insulated. Inputs with duplicate or
conflicting surface assignments, missing materials/capacities, unknown tags,
nonpositive capacities or nonfinite values are rejected.

## Verification

The `nova_gmsh_transient_heat_3d_verify` CTest executes a two-material
insulated wall with spatially uniform volumetric capacity `C=2` and source
`Q=4`: `T(t)=300+2t` is exact, irrespective of the differing conductivities.
The test reads final temperatures **from the exported VTK file** at `t=0.5`,
checks the discrete energy residual, and tests anisotropic materials with
Dirichlet/convection conditions and malformed inputs.

Verification is limited to the documented fixed-coefficient P1 thermal model.
There is no adaptive 3D time stepping, time-varying boundary data, radiation,
nonlinear material behavior, or transient restart CLI yet.
