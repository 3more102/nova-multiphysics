# Commercial readiness assessment

Release decision for NOVA 1.0 Commercial Edition: **NO-GO**. The assessed workflow is a stationary, homogeneous thermal box. The following weights apply only to this Engineering Alpha workflow. Give the full weight only when every named condition is evidenced. Partial work earns zero. Each gate remains mandatory independently of any subtotal.

| Gate 1 criterion for the thermal box | Weight | Testable acceptance | Earned |
| --- | ---: | --- | ---: |
| Working model representation | 2 | Stable geometry, mesh, material and study IDs, checked references and physical inputs, selections surviving remeshing | 2 |
| Native persistence | 2 | Version check, deterministic round trip, malformed-input rejection and preserved previous file on failure | 2 |
| Geometry to mesh | 2 | Parametric box, conforming tetrahedra, six selectable sides and successful remesh checks on all axes | 2 |
| Solver execution | 3 | Installed standalone process, exact wall temperature and heat-flow checks, refinement and heat-balance thresholds | 3 |
| Desktop authoring | 2 | Graphical create/edit/save/reopen and run of the supported study | 2 |
| Interactive result inspection | 1 | A user can select a location and inspect a computed field inside the application | 1 |

Earned **12/12** for this workflow. This score describes Engineering Alpha evidence for the thermal box. Full-platform maturity and commercial release remain subject to broader scope and independent gates.

| Area | Evidence on this development line | Open acceptance work |
| --- | --- | --- |
| Solver | Sparse CSR, CG/Jacobi/IC(0), GMRES and BiCGSTAB; true-residual reporting in the thermal runner. Cavity eigenanalysis uses a bounded reference solve. | Sparse eigenanalysis and larger workload qualification |
| Physics | Steady and transient conduction, electrostatics, small-strain isotropic elasticity, Maxwell/cavity and optical FEM have regression suites. Native file workflow covers steady homogeneous heat. | Intended-use boundaries, comparisons with independent solvers and experimental validation |
| Geometry | One parametric box with persistent side selections in a native project | General solid topology, repair and STEP interoperability |
| Meshing | Structured box tetrahedra and tagged Gmsh import in the SDK | Geometry-driven meshing of general solids and broader quality qualification |
| GUI and results | Installed Qt 6 editor, separate solver process, cancel, three orthogonal nodal planes, interactive kelvin probe and VTK export | Arbitrary-mesh 3D viewport, project history, usability evaluation |
| Multiphysics | Dedicated DC/Joule/thermal SDK workflows, including nonlinear and transient variants | File-driven coupled model and independent intended-use validation |
| Verification | Analytical tests, conservation, convergence, negative tests and installed consumers. The heated cube has a bounded continuum mean, an independent cell-centred finite-volume comparison, a full-field nodal implementation check, inversion-symmetry checks and a center-peak location check. | External engineering evaluation, experimental evidence and a continuum error bound for local peak temperature |
| Security and reliability | Size-limited project grammar, atomic file publication, CI sanitizer job | Distribution SBOM, dependency review, parser robustness campaign and signed-update verification |
| Licensing | Repository documents pending project-wide terms. Qt for Python is an optional desktop dependency. Core project I/O uses the C++ standard library. | Owner and legal approval for distribution, provenance and dependency audit |
| Distribution | Installable CMake SDK, runner, sample, optional desktop wheel, relocated shared-library runner and installed-consumer tests | Clean-machine installer, upgrade/uninstall tests and release credentials |
| Documentation | Format contract, equations, independent checker, benchmark and SDK installation guide | Product manuals, support and release maintenance procedures |

See the [pinned wall verification report](verification/thermal-project-v1.md), [heated cube comparison](verification/heated-cube.md) and [desktop workflow tests](THERMAL_DESKTOP.md). Gate 2 needs clean installation and external evaluation. Gates 3 and 4 additionally require platform qualification, completed security and licensing reviews, acceptance testing and owner approval for distribution.

## Engineering checkpoints, 25 November 2026

| Subsystem | Evidence submitted for review | Remaining release gate |
| --- | --- | --- |
| Adaptive tetrahedral diffusion | [PR #252](https://github.com/3more102/nova-multiphysics/pull/252) fixes an initially satisfied tolerance, reuses intermediate solves, and incorporates the two-generation mesh checks from #163. [Repository and release-gate inventory](https://github.com/3more102/nova-multiphysics/blob/8a498496b1a3418c13be3294d221ae56506b67bc/docs/verification/release-gates-2026-11-25.md) | Independent review and intended-use verification. The residual is a heuristic indicator. |
| Native tensor materials | [PR #249](https://github.com/3more102/nova-multiphysics/pull/249) checks all six components through saved projects with three orthogonal temperature gradients. A disposable combination of the two PRs passed 43 Linux CTest entries and six Qt workflows. | Review of the schema and distribution compatibility |
| Stationary thermo-optics | A [free-temperature equilibrium](THERMO_OPTIC_3D.md#coupled-equilibrium-with-four-free-thermal-nodes) compares six complex Maxwell edge integrals, four unconstrained temperatures, absorption and convection to an independent scalar solution. | Spatial convergence for coupled fields and intended-use validation |

The source version remains 0.8.0. No tested end-user installer, signed release or approved redistribution terms have been established by these checkpoints.
