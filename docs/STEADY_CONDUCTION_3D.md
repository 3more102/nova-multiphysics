# Steady conduction through tagged 3D solids

## Reproduce the simulation

From the repository root on a single-configuration CMake build:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
./build/nova_gmsh_diffusion_3d examples/gmsh_diffusion_3d/two_layer_wall.msh --material 11 2 --material 22 6 --temperature 101 400 --convection 102 12 300 --vtk wall.vtk
python verification/verify_tagged_conduction_3d.py ./build/nova_gmsh_diffusion_3d examples/gmsh_diffusion_3d/two_layer_wall.msh
```

On a Visual Studio multi-configuration build the executable is
`build/Release/nova_gmsh_diffusion_3d.exe`. Open `wall.vtk` as a legacy
unstructured grid. It has nodal `temperature`, cell `material_id` and cell
`heat_flux` vectors in global Cartesian coordinates. The verifier uses the
Python standard library. CTest registers it when a Python interpreter is available.
The C++ library, CLI, unit tests and wall smoke test run independently of Python.

## Model input

| Argument | Physical meaning in SI units |
| --- | --- |
| `--material volume_id k` | Isotropic conductivity in W/(m·K). Supply each volume ID once, including 0 for untagged tetrahedra when present. |
| `--material-tensor volume_id xx yy zz xy xz yz` | Symmetric positive-definite 3D conductivity tensor in global coordinates, W/(m·K). Six independent components; use either scalar or tensor for each volume ID. |
| `--temperature face_id T` | Prescribed nodal temperature in K on that surface. |
| `--flux face_id q_out` | Outward heat flux in W/m². Positive removes heat. |
| `--convection face_id h T_inf` | Outward flux equals `h (T − T_inf)`, with `h ≥ 0` in W/(m²·K). |
| `--source Q` | Uniform volumetric heating in W/m³. Default zero. |
| `--vtk output.vtk` | Optional file containing the solved nodal and cell fields. |

Use meters for mesh coordinates. Supply conforming four-node tetrahedra with shared
interface vertices. Surface triangle records identify *exterior* faces. An interior
physical triangle record is rejected. Omit interior surface records and identify
the adjoining solids by their volume IDs. Unassigned exterior faces are insulated.
The executable rejects absent, duplicate and conflicting assignments. The shorthand
`nova_gmsh_diffusion_3d mesh.msh tag` runs the unit-conductivity, unit-source example
with zero temperature on the selected tag.

## Discretization

For each tetrahedron, the diffusion matrix is the integral of
`grad(N_i) · K grad(N_j)`. Four positive degree-two volume points integrate
conductivity and source. Three degree-two face points integrate heat flux and
convection. Convection adds an `h N_i N_j` face matrix and an `h T_inf N_i` load.
A prescribed outward flux contributes `−q_out N_i` to the right-hand side.
Symmetric elimination lifts prescribed temperatures into the free equations.
NOVA assembles CSR and solves with Jacobi-preconditioned CG. Tagged mode uses a
relative true-residual tolerance of 1e-12 with zero absolute tolerance. An
unanchored connected component needs balanced loads and receives a volume-mean-zero
reference temperature.

The C++ problem accepts a spatial scalar callback, a spatial symmetric
positive-definite tensor callback, or `ScalarDiffusion3DProblem::cell_diffusivity`.
The vector follows `mesh.elements` ordering. An alternative
`ScalarDiffusion3DProblem::cell_diffusivity_tensors` holds one
`SymmetricTensor3D` per tetrahedron. The two cell arrays and either spatial
material callback are mutually exclusive. Both the linear solver and flux
postprocessor validate tensor positive definiteness. For mixed CLI material
types, scalar values are promoted to diagonal tensors, preserving the physical
isotropic conductivity.

Run a mixed isotropic/anisotropic tagged wall (here off-diagonals are zero so
the one-dimensional analytical reference remains valid):

```sh
./build/nova_gmsh_diffusion_3d examples/gmsh_diffusion_3d/two_layer_wall.msh \
  --material-tensor 11 2 9 4 0 0 0 --material 22 6 \
  --temperature 101 400 --convection 102 12 300 --vtk anisotropic-wall.vtk
```

General nonzero cross terms produce transverse heat flux under an axial
temperature gradient. Exterior insulation and the prescribed boundaries must
therefore be physically consistent; no one-dimensional solution is assumed in
that case. Assembly and flux recovery select the
same material on each tetrahedron. Continuous P1 interpolation enforces temperature
continuity across a bonded interface. The recovered flux is the cell average of
`−K grad(T_h)`. Its tangential components can jump across materials. General raw
P1 face fluxes need a separate equilibration procedure for local conservation.

## Analytical verification

The supplied 1 m wall has unit cross section and two 0.5 m layers. The hot surface
is 400 K. The cold surface exchanges heat with 300 K ambient at 12 W/(m²·K).
The series resistance per unit area is `0.5/2 + 0.5/6 + 1/12`. This gives an
axial flux of 240 W/m², an interface temperature of 340 K and an outer surface
temperature of 320 K. Swapping the conductivities moves the interface to 380 K.
Reference: [COMSOL Heat Transfer Module User’s Guide, thermal resistances in series](https://doc.comsol.com/6.3/doc/com.comsol.help.heat/HeatTransferModuleUsersGuide.pdf).

The following GCC 14.2 Release results were read from the 12-tetrahedron VTK file.
Temperature and flux gates allow margin above the 1e-12 linear stopping tolerance
for gradient recovery and cross-platform floating-point differences. For a zero
reference, the absolute residual is the acceptance metric.

| Quantity | Reference | Computed | Absolute error | Relative error | Absolute limit |
| --- | ---: | ---: | ---: | ---: | ---: |
| Worst interface temperature, K | 340 | 340.00000000000011 | 1.14e-13 | 3.34e-16 | 2e-8 |
| Worst cold-surface temperature, K | 320 | 320.00000000000011 | 1.14e-13 | 3.55e-16 | 2e-8 |
| Worst axial heat flux, W/m² | 240 | 239.99999999999864 | 1.36e-12 | 5.68e-15 | 2e-7 |
| Hot-face outward rate, W | −240 | −240.00000000000023 | 2.27e-13 | 9.47e-16 | 2e-7 |
| Cold-face outward rate, W | 240 | 240.00000000000068 | 6.82e-13 | 2.84e-15 | 2e-7 |
| Net outward heat rate, W | 0 | 4.55e-13 | 4.55e-13 | — | 2e-7 |

A second run sets a uniform source of 120 W/m³ and outward flux of 240 W/m² at
the right boundary. The analytical flux is `q_x(x) = 120 + 120 x`. The temperature
starts at 400 K with slope `−q_x/k` in each layer. An independent positive
Duffy-transformed tensor Gauss rule integrates the squared temperature error
(degree four) and the squared vector flux error over every exported tetrahedron.

| Tetrahedra | Temperature L2 error | Order | Flux L2 error | Order |
| ---: | ---: | ---: | ---: | ---: |
| 48 | 1.007988826 | — | 16.86907852 | — |
| 384 | 0.259021484 | 1.9603 | 8.560864074 | 0.9786 |
| 3,072 | 0.065282435 | 1.9883 | 4.312253592 | 0.9893 |

The unit-volume domain makes the L2 errors numerically equal to RMS errors in K
and W/m². Acceptance requires orders of at least 1.7 and 0.85 on both refinements
and fine-grid errors below 0.08 K and 5 W/m². Both the analytical and the
refinement tables passed. The end-to-end verifier also checks a mixed tensor/scalar tagged wall against
the same analytical temperature, axial flux, transverse flux and heat-rate
targets, and confirms invalid material tensors cannot replace a prior VTK
result. Separate C++ regressions exercise per-cell non-diagonal tensor
stiffness with winding reversal, two anisotropic layers with exact affine
solutions, continuous normal flux and discontinuous transverse flux, and
rejection of non-SPD or ambiguous material data.

Separate C++ regressions exercise an oblique interface
with conserved normal flux and discontinuous tangential flux, symmetry and the
constant mode of each element matrix, reversed element orientation, bad material
data and incompatible pure-Neumann loads. The CLI verifier also checks that failed
runs preserve an existing result.

## Performance measurements

Append `--benchmark` to the Python command above for a three-run median of the
complete imported-mesh solve and VTK export at each size. Numbers below use GCC
14.2 Release on an x86-64 Linux AMD EPYC 9V74 host. Wall times below are medians
of five measured runs per size. Peak resident memory comes from per-process
`wait4`. CG times are medians of three separate instrumented runs. File I/O is
included in wall time. The solver runs serially.

| Tetrahedra | DOFs | Constrained CSR nonzeros | CG iterations | CG time, ms | End to end, s | Peak RSS, MiB |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3,072 | 729 | 4,295 | 67 | 0.56 | 0.0204 | 4.27 |
| 24,576 | 4,913 | 31,503 | 128 | 5.12 | 0.1319 | 20.34 |
| 82,944 | 15,625 | 103,127 | 183 | 22.91 | 0.5300 | 63.71 |

On the largest mesh the FEM solve, including geometry checks, assembly, CSR
construction and CG, took a median 0.213 s. CSR construction took 0.026 s.
The constrained triplet buffer had 1,319,665 entries before CSR coalescence.
Profiling identified five complete mesh topology traversals. The importer now
obtains validation and exterior faces in the same traversal, reducing the full
workflow to four. In alternating before/after measurements, the largest case moved
from 0.5655 s to 0.5300 s. The smallest case varied from 0.0179 s to 0.0204 s.
These timings characterize this mesh family and host.

Each P1 element contributes 16 stiffness triplets. CSR matrix-vector work scales
with the stored nonzeros per iteration. The current triplet buffer, ordered maps
and geometry validation account for memory beyond the final CSR arrays. Element
integration is a candidate for parallel execution if deterministic shared-node
assembly is addressed. The present verification covers steady bonded conduction,
conforming linear tetrahedra and the stated mesh sizes. Thermal contact resistance,
experimental validation and industrial-scale performance remain open.
