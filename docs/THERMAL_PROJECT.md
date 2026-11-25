# Native thermal project, schemas 1 and 2

The installed C++20 `nova-run` process reads a self-contained SI model. It builds the mesh, assigns the material and boundary conditions, solves stationary heat conduction, and returns temperatures and heat rates. The model and codec are also available to SDK clients in `nova/model/thermal_project.hpp` and `nova/io/thermal_project.hpp`.

## Run the reference wall

~~~sh
nova-run check convecting_wall.nova
nova-run show convecting_wall.nova
nova-run save-as convecting_wall.nova saved.nova
nova-run solve saved.nova --vtk wall.vtk
~~~

The `show` operation emits the canonical versioned model to standard output. [Desktop authoring and field probes](THERMAL_DESKTOP.md) · [Numerical verification report](verification/thermal-project-v1.md).

The SDK installs the runner in `bin` and the project examples, including `anisotropic_wall.nova`, in `share/nova/projects` under the default install layout. A CMake consumer can use `NOVA::Runner` and `NOVA_PROJECT_EXAMPLE_DIR`. Open the VTK unstructured grid in ParaView to inspect point temperature in K and cell heat flux in W/m². Material tag 1 designates the single homogeneous solid.

## Model contract

Use printable ASCII, spaces or tabs between fields, and LF or CRLF line endings. Blank lines and comments beginning with `#` are accepted. Each non-comment line has the indicated fields in order. Example input:

~~~text
NOVA_PROJECT 1
project convecting-wall
geometry wall box 0.5 0.5 0.5
mesh tets wall structured_tet 8 3 2
material homogeneous isotropic 2
study steady steady_heat tets homogeneous
source steady 0
solver steady cg jacobi 2000 1e-13 1e-14
boundary steady wall xmin temperature 400
boundary steady wall xmax convection 12 300
end
~~~

The distributed `examples/projects/anisotropic_wall.nova` uses a synthetic material for an exact field test. A schema-2 anisotropic project uses the same model records with the tensor material form:

~~~text
NOVA_PROJECT 2
project anisotropic-wall
geometry wall box 0.5 0.25 0.125
mesh tets wall structured_tet 6 2 1
material benchmark-solid anisotropic 2 3 4 0.25 -0.15 0.10
study steady steady_heat tets benchmark-solid
source steady 0
solver steady cg jacobi 2000 1e-13 1e-14
boundary steady wall xmin temperature 400
boundary steady wall xmax convection 12 300
boundary steady wall ymin outward_flux -37.5
boundary steady wall ymax outward_flux 37.5
boundary steady wall zmin outward_flux 22.5
boundary steady wall zmax outward_flux -22.5
end
~~~

| Record | Meaning and units |
| --- | --- |
| project | Unique persistent project ID |
| geometry | Geometry ID, box, three positive axis lengths in m, with the minimum corner at the origin |
| mesh | Mesh ID, geometry ID, structured tetrahedra, three positive subdivision counts |
| material | Material ID and either positive isotropic conductivity or a symmetric positive-definite tensor in W/(m·K) |
| study | Study ID, stationary heat conduction, mesh ID and material ID |
| source | Study ID and spatially uniform signed volumetric generation in W/m³ |
| solver | Study ID, CG, preconditioner (`none`, `jacobi`, `ic0`), iteration budget, relative and absolute residual tolerances |
| boundary | Study ID, geometry ID, one of `xmin xmax ymin ymax zmin zmax`, then `temperature T` (K), `outward_flux g` (W/m²) or `convection h T_ambient` (W/(m²·K), K) |
| end | Required final record |

Schema 1 retains the original isotropic material record and is still written for isotropic projects. Schema 2 adds `material ID anisotropic Kxx Kyy Kzz Kxy Kxz Kyz`, where the six finite entries define a symmetric conductivity tensor in the global Cartesian frame, in W/(m·K). The tensor must be numerically positive definite under the same scaled Cholesky check used by FEM assembly. Its full tensor is used for both the weak form and Fourier heat flux; off-diagonal terms are not discarded. Schema 2 also accepts isotropic material records. The box, mesh, study, material and project use distinct identifiers of up to 64 ASCII characters, beginning with a letter and continuing with letters, digits, underscores or hyphens. Boundary records can appear in any order after the solver. Each side accepts one assignment. Unassigned sides are insulated. At least one fixed temperature or positive convection coefficient anchors absolute temperature. Intersecting fixed sides must agree. Changing the dimensions or the mesh density preserves side selections. Parameters are model inputs, and the example uses constructed verification values.

`check` parses the record, validates physical inputs and references, and constructs and validates the mesh. `save-as` writes canonical ASCII with round-trip double precision, reopens the staged record and replaces the destination after successful read-back. The format rejects unknown versions, missing or extra fields, trailing data, invalid references, and nonfinite numbers. Project streams are limited to 64 KiB, generated meshes to 196,608 tetrahedra and CG to 100,000 iterations. These are file-runner resource limits. The underlying numerical SDK has its own interfaces. Relative tolerance must be in [0, 1), absolute tolerance nonnegative, and at least one must be positive.

A successful solve emits a `NOVA_RESULT 1` header followed by named scalar fields, and optionally a VTK result. Positive heat rates are outward. `residual` is the CG residual of the constrained linear system. `max_free_residual_W` and `balance_defect_W` come from reassembly of the thermal weak form. Failed solves return a nonzero process status and leave any previous VTK result intact. Project and VTK writes use a same-filesystem staging directory. Recovery from whole-machine power loss or an untrusted shared output directory is outside this contract.

## Formulation and verification

The runner uses NOVA's existing P1 tetrahedral diffusion assembly for `−div(K grad T) = Q` and Fourier heat flux `q = −K grad T`. Schema 1 uses `K = kI`; schema 2 can use a homogeneous symmetric positive-definite anisotropic tensor. Exterior Neumann data prescribe `q·n = g`; convection prescribes `q·n = h(T−T∞)`. The balance is generated power minus prescribed, convective and fixed-temperature outward power. The adapter calls the FEM routines for flux recovery and reaction accounting.

For the supplied wall, `q_x = (T_hot−T∞)/(L/k + 1/h)`, giving 300 W/m², 75 W through the 0.25 m² cross section, and a 325 K cold surface. The independent verifier reads every exported nodal temperature and cell flux. For the heated wall it uses `T(x) = T_hot + ((QL−g)x − Qx²/2)/k` and `q_x(x) = g − Q(L−x)`. It evaluates cell-centroid temperature and vector-flux RMS errors on four uniform refinements. These sampled norms use equal-volume tetrahedra. The test checks boundary reactions and power separately. [Plane-wall heat equation derivation, Ankara University](https://acikders.ankara.edu.tr/pluginfile.php/206354/mod_resource/content/0/ENE%20305335-Chapter%203%20continues-1.pdf).

~~~sh
ctest --test-dir build -R 'nova_thermal_project|nova_unit_tests' --output-on-failure
python verification/verify_thermal_project.py generated build/nova-run examples/projects/convecting_wall.nova
python verification/benchmark_thermal_project.py build/nova-run examples/projects/convecting_wall.nova 8 16 32
~~~

The benchmark records end-to-end solve time without VTK export, process peak RSS (Linux KiB), CG iterations and true residual, cold-surface and power error, and energy defect. Each row starts a new process. Supported geometry is a single axis-aligned box with one homogeneous material. Schema 1 is isotropic; schema 2 adds one constant anisotropic tensor. Transient studies, external meshes, spatially varying materials, electric coupling and material catalogues are not yet part of the native project format. The Qt desktop editor can author, save, reopen and solve both schema-1 isotropic and schema-2 constant-anisotropic projects. Legacy SDK APIs and checkpoint formats keep their own version contracts.

The [internally heated cube verification](verification/heated-cube.md) checks genuinely three-dimensional conduction against a bounded continuum series, a nodal SciPy implementation check and an independent cell-centred finite-volume calculation. Open `internally_heated_cube.nova` from the installed projects directory to reproduce the study.


### Verification of all six tensor entries

The `anisotropic` verification suite solves the distributed x-directed wall
and generates two additional projects from the same material and mesh.
The exact temperatures are 400 − 150x, 400 − 150y and 400 − 150z K, with
coordinates in metres. Each selects a different column of the fixed tensor.
The prescribed lateral fluxes and convective ambient enforce the same affine
field throughout the box.

| Temperature gradient direction | Exact heat flux (W/m²) | Convective removal (W) |
| --- | --- | ---: |
| x | (300, 37.5, −22.5) | 9.375 |
| y | (37.5, 450, 15) | 28.125 |
| z | (−22.5, 15, 600) | 75 |

For each solve the checker evaluates all 42 nodal temperatures, all 72 cell
fluxes, the opposing thermal reactions and net energy balance. The absolute
thresholds are 2e-8 K for temperatures, 2e-8 W/m² for vector flux error, and
2e-8 W for integrated powers. These manufactured cases exercise serialization,
assembly and postprocessing for a homogeneous tensor on a box.
