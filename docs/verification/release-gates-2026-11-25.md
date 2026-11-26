# NOVA 1.0 release-gate inventory, 25 November 2026

**Commercial release decision: blocked.** CMake declares version 0.8.0. The audited main commit is [ed41d9a](https://github.com/3more102/nova-multiphysics/commit/ed41d9afa6794fe8c7cf0a714a00bed55bc8a4e4). Results described for PR #252 belong to its integration branch until it merges. The release gate for each row below remains incomplete. A specific benchmark can pass while qualification for a subsystem remains open.

## Backlog and build evidence

The opening snapshot had 262 remote branch heads, 74 open pull requests and zero open issues. A three-way Git merge-tree check against the pinned main found 30 clean merges and 44 with conflicts. None of the 30 clean merges produced the main tree or the same tree as another clean merge. The check uses source trees only. A clean merge still needs semantic review and CI on the resulting tree.

| Merge-tree result | PR numbers at the opening snapshot |
| --- | --- |
| Clean | 252, 249, 248, 242, 239, 238, 237, 234, 233, 223, 221, 211, 209, 208, 206, 205, 200, 198, 176, 173, 169, 168, 167, 166, 165, 163, 162, 160, 158, 156 |
| Conflicts | 241, 240, 236, 235, 232, 229, 228, 227, 226, 224, 222, 220, 219, 218, 217, 216, 215, 214, 213, 210, 207, 204, 203, 201, 199, 192, 190, 188, 187, 171, 164, 161, 159, 155, 131, 17, 14, 13, 12, 10, 9, 7, 6, 5 |

- [PR #252](https://github.com/3more102/nova-multiphysics/pull/252) now incorporates the two-generation mesh transfers, boundary winding and material-volume checks from [#163](https://github.com/3more102/nova-multiphysics/pull/163). Its installed SDK smoke test exercises the bounded-driver use case from [#160](https://github.com/3more102/nova-multiphysics/pull/160). The latter proposes a different public history type and needs an API decision before closure.
- [Issue #253](https://github.com/3more102/nova-multiphysics/issues/253) records premature refinement of an already-accepted mesh. Its regression failed on the original #252 head. The correction and the regression are in the same PR.
- [PR #249](https://github.com/3more102/nova-multiphysics/pull/249) changes the project schema, editor and runner for anisotropy. It is a separate file-format and product compatibility decision from the solver change.
- Overlapping candidates include zero-permittivity tests (#237/#238), oblique multilayers (#192/#216/#222), slab modes (#198/#215), rings (#213/#214/#217/#219) and dispersive materials (#199/#200/#201/#206). Preserve distinct physical cases, interfaces and tests when reconciling these branches.
- The baseline [main CI run](https://github.com/3more102/nova-multiphysics/actions/runs/45470001906) passed. Subsequent branch results are available in [PR checks](https://github.com/3more102/nova-multiphysics/pull/252/checks). No line/branch coverage measurement was established in this checkpoint.

Reproduce triage with `git fetch --prune origin` in a disposable clone, GitHub's paginated open PR list, and `git merge-tree --write-tree <base-sha> <pr-head-sha>` for each pair. The check identified integration conflicts without modifying any development branch.

## Platform inventory

| Release gate | Classification | Implemented evidence | Outstanding release evidence or missing capability |
| --- | --- | --- | --- |
| Architecture | Partially implemented | C++20 core, installed targets, separate desktop/runner, native schema and module headers | Stable API, migration commitments and a unified multiphysics model |
| FEM and algebra | Implemented, insufficiently verified for general use | 1D P1, 2D P1/P2 triangles, 3D P1 tetrahedra, CSR, CG/Jacobi/IC(0), GMRES, BiCGSTAB, Newton, and local tetrahedral refinement | Quad/hex formulations, mesh coarsening, general sparse direct solve and larger independent qualification suite |
| Thermal | Implemented, insufficiently verified for general use | Mixed-BC stationary and transient conduction, tensors, energy audits, nonlinear electrothermal material laws, wall and heated-cube checks | Radiation, contact resistance, qualified temperature-dependent material ranges and intended-use validation |
| Electromagnetics | Experimental | Electrostatics, terminal models, driven Nédélec tetrahedra, energy/power audits and bounded PEC cavity reference | General magnetostatic/eddy-current workflows, ports, open boundaries and scalable eigensolve |
| Photonics | Experimental | 1D complex Helmholtz FEM, normal/oblique optical checks and planar analytical references | General vector waveguide modes, resonator studies, device-level scattering and measured comparisons |
| Multiphysics and structures | Experimental | Joule heating, transient electrothermal, Maxwell-to-thermal absorption, stationary self-consistent thermo-optics, small-strain linear elasticity | Reusable general coupling framework, qualified coupling to structural mechanics and independent application-level validation |
| Geometry | Partially implemented | Native box and named sides | STEP/IGES/B-rep import, repair, topology identity and geometric operations. Open CASCADE is a candidate C++ modeling/data-exchange kernel; production integration needs a versioned adapter and representative CAD corpus. [Vendor overview](https://occt3d.com/open-cascade-technology/index.html) |
| Meshing | Partially implemented | Structured box tetrahedra, triangular meshes, tagged Gmsh ASCII 2.2/4.1 inputs, mesh-quality checks, refinement and transfer | Automatic volume mesh generation from general solids, curvature sizing and boundary layers |
| GUI and project workflow | Partially implemented | Qt thermal editor, model tree, save/reopen, solver subprocess, cancellation, field probes and VTK export | General model authoring, undo/redo, desktop autosave/recovery and arbitrary-mesh visualization. Solver restart files cover narrower cases. |
| Materials | Partially implemented; catalog missing | User-entered model coefficients, spatial callbacks, scalar/tensor cell data and optical reference laws | Versioned, source-attributed catalog with validity ranges and interoperable materials import/export |
| Performance and HPC | Partially implemented | Sparse C++ core and Linux wall runtime/RSS harness | Repeated scaling with controlled hardware, large-model memory limits, threaded and distributed qualification |
| Reliability and V&V | Partially implemented | Analytical cases, convergence tests, SDK/CLI and desktop tests, checkpoints and sanitizer CI | Intended-use qualification, independent review, failure/recovery coverage and external validation |
| Security and supply chain | Partially implemented | Input validation, bounded project parsing, sanitizer checks, read-only workflow permissions | Shipped-binary SBOM, dependency findings/disposition, integrity signatures and update verification |
| Installation | Partially implemented | CMake install/export, relocatable shared SDK, installed consumers on CI and an optional Python desktop package | Clean Windows 11 and Linux end-user installer, upgrade, uninstall and recovery qualification |
| Licensing | Blocked | README records pending project-wide terms. No top-level project license found at the audited revision. | Owner-approved terms, contribution provenance, component notices, legal review and a defined redistribution configuration |
| Documentation | Partially implemented | SDK, thermal project/desktop, FEM, solver and benchmark guides in docs/ | Cross-checked full manuals, complete API reference, user support and release operations |
| Versioned release | Blocked | Source version and CI for development packages | Signed reproducible installer, complete release evidence and approval |

### Dependency evidence

CMake links the core without a third-party solver target. The optional desktop declares PySide6-Essentials in its Python manifest. The CI-tested Linux wheels for PySide6-Essentials and shiboken6 were both 6.8.3. Qt documents community and commercial distributions and separate third-party components. Record the exact platform wheels, Qt libraries, plugins, Python/compiler runtimes and notices for each planned installer before making a redistribution decision. SciPy 1.17.0 is used by the independent verification job. The CI workflow also installs system graphics libraries and references checkout/setup-python by major version tags. The transitive license and vulnerability inventory of a final installer remains open. [Qt for Python licensing](https://doc.qt.io/qtforpython-6/commercial/index.html) · [Qt third-party components](https://doc.qt.io/qtforpython-6/licenses.html) · [Desktop manifest](../../desktop/pyproject.toml) · [CI workflow](../../.github/workflows/ci.yml).

### Measured engineering checkpoint

Debian 13 x86-64, GCC 14.2 Release. The revised source passed 42 CTest entries, including the main test executable (414 C++ cases). Four installed SDK checks passed. Five offscreen Qt workflows passed. These counts cover different suites and should be read separately.

| Numerical case | Measurement |
| --- | --- |
| Manufactured adaptive wall | 48 to 314 tetrahedra; exact energy-seminorm error 0.408248290464 to 0.371639993999. Minimum mesh mean ratio moved from 0.75595 to 0.29075. |
| Layered anisotropic wall with Robin and signed Neumann data | Independent series-resistance solution, axial flux 8 and transverse fluxes 2 and −1. Checked two aligned grids and two algebraic tolerances. See the [acceptance limits](../ADAPTIVE_DIFFUSION_3D.md). |
| Installed heated cube | At 16 subdivisions, the independent SciPy nodal implementation differed by at most 8.6402e-12 K over 4,913 nodes. The *continuum* upper bound on mean-temperature error was 0.013858456 K. The separate cell-centered finite-volume reference achieved an upper bound of 0.00094203778 K at 64 subdivisions. [Problem and reference definitions](heated-cube.md) |

Single-process wall benchmark on an AMD EPYC 9V74 virtualized Linux host, one observation per size. Elapsed time includes runner startup and output. These are baseline observations without isolated-hardware or repeat-run uncertainty estimates.

| Subdivisions per axis | Tetrahedra | Wall seconds | Peak RSS (KiB) | CG iterations | Heat balance defect (W) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 3,072 | 0.016163 | 10,624 | 71 | 3.83e-12 |
| 16 | 24,576 | 0.108199 | 19,092 | 137 | 2.24e-11 |
| 32 | 196,608 | 1.364238 | 137,244 | 253 | 1.60e-10 |

From the source root after configuring a Release build and installing to `sdk`:

~~~sh
ctest --test-dir build --output-on-failure
cmake -S tests/package_consumer -B consumer -DCMAKE_PREFIX_PATH="$PWD/sdk"
cmake --build consumer
ctest --test-dir consumer --output-on-failure
build/nova_adaptive_diffusion_3d
python verification/verify_thermal_project.py cube sdk/bin/nova-run sdk/share/nova/projects/internally_heated_cube.nova --scipy
python verification/benchmark_thermal_project.py sdk/bin/nova-run sdk/share/nova/projects/convecting_wall.nova 8 16 32
~~~
