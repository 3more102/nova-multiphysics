# Conservative face-flow projection for 3D steady conduction

The C++20 API `nova::fem::reconstruct_conservative_heat_flow_3d` gives every
tetrahedral face exactly one signed integrated heat rate, positive outward
from its owner cell. Both incident tetrahedra use the same interior rate with
opposite signs. The source is integrated with the same positive degree-two
tetrahedral rule as the steady FEM solver.

## Reconstruction

1. Evaluate the existing P1 conductivity-dependent cell-average fluxes.
2. Use the mean of both incident raw cell-flux projections on every interior
   face. For prescribed exterior Neumann faces, integrate their given outward
   flux; for Robin faces, integrate `h (T - T_ambient)` with the face
   degree-two rule. Insulated faces remain exactly zero.
3. Calculate each cell defect `source_heat - sum(outward_face_heat)`.
4. For each connected tetrahedral-cell component, transfer deficits along
   a deterministic spanning tree of interior faces. A Dirichlet exterior
   face, if present, absorbs the final component deficit.
5. With no Dirichlet face on a component, source and prescribed exterior
   outflow must be globally compatible (within integration-roundoff
   tolerance); otherwise the postprocessor rejects the configuration.

This is an **algebraic conservative projection**, not an H(div)-conforming
RT0 finite-element reconstruction or minimum-norm equilibration. Its
corrections and especially the split of reaction heat among Dirichlet faces
are not unique. The tree can concentrate corrections on individual faces.
The result is a conservative face-flow *diagnostic*; it does not make a
nonconverged temperature field a correct physical solution. Use
`evaluate_diffusion_energy_balance_3d` to audit the assembled FEM free-node
residuals independently and examine `maximum_absolute_face_correction`
as a warning of mismatch between raw P1 fluxes and the conservative rates.

Explicit interior Dirichlet **vertex** constraints do not create a physical
exterior outlet. Where these induce nonzero reaction heat but no Dirichlet
boundary face is declared, the reconstruction rejects incompatible loads
rather than inventing a surface reaction.

## CLI

With a Gmsh 2.2 mesh and the existing material/boundary options:

```sh
nova_gmsh_diffusion_3d two_layer_wall.msh \
  --material 11 2 --material 22 6 \
  --temperature 101 400 --convection 102 12 300 \
  --balance-report --local-balance-report
```

The CLI prints the maximum raw and reconstructed cell heat defects,
largest face correction, generated heat and total exterior heat outflow.
The per-face rates and per-cell arrays are available through the public
C++ SDK API. Existing default solver output is unchanged.

## Verification

`tests/test_conservative_face_flow_3d.cpp` covers analytic single-tetra
source/Neumann conservation, an interior interface with distinct anisotropic
materials, a sourced multi-material Robin/Dirichlet wall, zero insulated
flow, and rejection of incompatible or malformed boundary conditions.
CMake registers the core test suite and a separate CLI smoke test.
