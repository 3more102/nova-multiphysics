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
`heat_flux` vectors in global Cartesian coordinates. With `--heat-source`
the VTK result additionally contains `heat_generation` in W/m³ per cell,
including the uniform `--source` contribution if specified. The verifier uses the
Python standard library. CTest registers it when a Python interpreter is available.
The C++ library, CLI, unit tests and wall smoke test run independently of Python.

## Model input

| Argument | Physical meaning in SI units |
| --- | --- |
| `--material volume_id k` | Isotropic conductivity in W/(m·K). Supply each volume ID once, including 0 for untagged tetrahedra when present. |
| `--temperature face_id T` | Prescribed nodal temperature in K on that surface. |
| `--flux face_id q_out` | Outward heat flux in W/m². Positive removes heat. |
| `--convection face_id h T_inf` | Outward flux equals `h (T − T_inf)`, with `h ≥ 0` in W/(m²·K). |
| `--source Q` | Uniform volumetric heating in W/m³. Default zero. |
| `--heat-source volume_id Q` | Signed volumetric source density in W/m³ for a selected volume ID. Can repeat for different IDs; unlisted IDs add zero. This density adds to `--source`. |
| `--vtk output.vtk` | Optional file containing the solved nodal and cell fields. |
| `--flux-report` | Print per-physical-surface triangle count, area (m²), and outward conductive heat rate (W), including untagged ID 0. |
| `--balance-report` | Print consistent FEM-integrated generation, physical Dirichlet reaction flow, Neumann and Robin heat flow, plus signed and relative conservation defect. |

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

The same physical-tag table can now supply different volumetric heat-generation
densities. For example, `--heat-source 11 120 --heat-source 22 360` heats the
two layers at 120 and 360 W/m³ respectively. Negative values represent cooling.
The optional `ScalarDiffusion3DProblem::cell_source` is additive to the
spatial source callback and follows exactly the tetrahedron element order.
Nonfinite source densities and incomplete per-cell source vectors are rejected.
Per-volume sources are applied at volume quadrature points, including in the
pure-Neumann compatibility test.

With the 400 K hot wall, 12 W/(m²·K) convection to 300 K and the
two conductivities 2 and 6 W/(m·K), the independent one-dimensional solution
has heat flux 144 W/m² at the hot wall, 204 W/m² at the material interface
and 384 W/m² at the cold surface; the interface temperature is 356.5 K
and the cold temperature is 332 K. The total volumetric generation is 240 W.
The Python verifier checks exported temperature/flux convergence on
successively refined meshes, rather than assuming raw P1 element flux
is pointwise exact where sources are present.

The C++ problem accepts a spatial scalar callback, a spatial symmetric
positive-definite tensor callback, or `ScalarDiffusion3DProblem::cell_diffusivity`.
The vector follows `mesh.elements` ordering. Assembly and flux recovery select the
same material on each tetrahedron. Continuous P1 interpolation enforces temperature
continuity across a bonded interface. The recovered flux is the cell average of
`−K grad(T_h)`. Its tangential components can jump across materials. General raw
P1 face fluxes need a separate equilibration procedure for local conservation.

## Tagged exterior boundary heat-flow diagnostic

Append `--flux-report` to the Gmsh example command to print signed outward
conductive heat rates for every surface physical ID. The report uses the *owner*
tetrahedron's cell-averaged raw P1 flux and the geometric outward face-area
vector. Face winding and tetrahedral node ordering do not affect the result.
The API `evaluate_tagged_boundary_flux_3d` also checks for duplicate faces,
incorrect vector lengths, nonfinite fluxes, and interior triangles passed
as exterior surfaces.

On the two-layer analytical wall, physical ID 101 (hot surface) has outward
rate **−240 W** and ID 102 (convection surface) has **+240 W**, each over
1 m². Untagged insulated faces have zero axial conduction; the net is zero
within numerical tolerance. Positive reported values mean heat exits the solid.

**Numerical limitation:** this diagnostic integrates the raw cell-average flux,
not the prescribed Neumann/convection load or the FEM Dirichlet reaction.
It is exact for constant conductivity and affine temperatures on each cell,
as in the analytical wall. With spatially varying conductivity, a nonzero
source, or general finite-element errors, raw P1 fluxes are not locally
equilibrated. Reported surface rates therefore are estimates and their sum
is **not** an enforced conservation identity. The volume-integrated source
and FEM reaction require a separate consistent balance computation before
the output can be used as a strict energy-conservation certificate.

## Consistent finite-element energy balance

Use `--balance-report` to request a physically signed global heat budget in
watts, independent of the optional raw `--flux-report`. Positive outward rates
remove heat from the solid. The diagnostic evaluates source density with the
**same tetrahedral quadrature used in assembly**, records prescribed Neumann
integrals, and computes Robin heat exchange from the originally sampled face
quadrature using the solved P1 face temperature. The Dirichlet contribution
comes from unconstrained nodal FEM reactions **before** rows are eliminated
or artificial pure-Neumann reference pins are applied. In the solver's weak
form, the physical Dirichlet outward heat rate equals minus `(A*T - b)`
summed over genuinely prescribed nodes.

The resulting discrete energy identity is

```text
Q_Dirichlet,out + Q_Neumann,out + Q_Robin,out - Q_generated = imbalance
relative_imbalance = abs(imbalance) / max(abs(each contribution))
```

An exactly zero reference uses relative imbalance zero. The reported quantities
are global integrals, **not** separately equilibrated local face fluxes.
The relative defect is a diagnostic of *discrete* global balance, not a
physical accuracy estimate: convergence of the temperature field, gradients,
and boundary fluxes still needs analytical or mesh-refinement verification.
The solve must also converge before interpreting the result. The optional
`ScalarDiffusion3DResult::energy_balance` is absent unless
`ScalarDiffusion3DProblem::compute_energy_balance` is enabled; normal
solves avoid storing additional uneliminated reaction rows.

For the source-free analytical wall, Dirichlet rate is −240 W, Robin rate
is +240 W, and the generated rate is zero. With
`--heat-source 11 120 --heat-source 22 360`, the total integrated source
is **240 W**. The conservation defect should remain small even when the
raw element-face flux diagnostic does not close due to discretization error.
The new tests exercise signed Dirichlet reactions, simultaneous spatial and
cellwise generation, prescribed Neumann flux, convection, and pure-Neumann
gauge exclusion.

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
refinement tables passed. Separate C++ regressions exercise an oblique interface
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
