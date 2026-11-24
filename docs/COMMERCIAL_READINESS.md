# Commercial readiness assessment

Release decision for NOVA 1.0 Commercial Edition: **NO-GO**. The assessed workflow is a stationary, homogeneous thermal box. The following weights apply only to this Engineering Alpha workflow. Give the full weight only when every named condition is evidenced. Partial work earns zero. Each gate remains mandatory independently of any subtotal.

| Gate 1 criterion for the thermal box | Weight | Testable acceptance | Earned |
| --- | ---: | --- | ---: |
| Working model representation | 2 | Stable geometry, mesh, material and study IDs, checked references and physical inputs, selections surviving remeshing | 2 |
| Native persistence | 2 | Version check, deterministic round trip, malformed-input rejection and preserved previous file on failure | 2 |
| Geometry to mesh | 2 | Parametric box, conforming tetrahedra, six selectable sides and successful remesh checks on all axes | 2 |
| Solver execution | 3 | Installed standalone process, exact wall temperature and heat-flow checks, refinement and heat-balance thresholds | 3 |
| Desktop authoring | 2 | Graphical create/edit/save/reopen and run of the supported study | 0 |
| Interactive result inspection | 1 | A user can select a location and inspect a computed field inside the application | 0 |

Earned **9/12** for this workflow. This score describes Engineering Alpha evidence for the thermal box. Full-platform maturity and commercial release remain subject to broader scope and independent gates.

| Area | Evidence on this development line | Open acceptance work |
| --- | --- | --- |
| Solver | Sparse CSR, CG/Jacobi/IC(0), GMRES and BiCGSTAB; true-residual reporting in the thermal runner. Cavity eigenanalysis uses a bounded reference solve. | Sparse eigenanalysis and larger workload qualification |
| Physics | Steady and transient conduction, electrostatics, small-strain isotropic elasticity, Maxwell/cavity and optical FEM have regression suites. Native file workflow covers steady homogeneous heat. | Intended-use boundaries, comparisons with independent solvers and experimental validation |
| Geometry | One parametric box with persistent side selections in a native project | General solid topology, repair and STEP interoperability |
| Meshing | Structured box tetrahedra and tagged Gmsh import in the SDK | Geometry-driven meshing of general solids and broader quality qualification |
| GUI and results | Versioned solver response and VTK temperature, heat flux and generation; external visualization | Desktop authoring, in-application field inspection and usability evaluation |
| Multiphysics | Dedicated DC/Joule/thermal SDK workflows, including nonlinear and transient variants | File-driven coupled model and independent intended-use validation |
| Verification | Analytical tests, conservation, convergence, negative tests and install-consumer suite | External engineering evaluation and experimental evidence |
| Security and reliability | Size-limited project grammar, atomic file publication, CI sanitizer job | Distribution SBOM, dependency review, parser robustness campaign and signed-update verification |
| Licensing | Repository documents pending project-wide terms. The core thermal project adds no third-party runtime dependency. | Owner and legal approval for distribution, provenance and dependency audit |
| Distribution | Installable CMake SDK, runner, sample and installed-consumer test | Clean-machine installer, upgrade/uninstall tests and release credentials |
| Documentation | Format contract, equations, independent checker, benchmark and SDK installation guide | Product manuals, support and release maintenance procedures |

Regression entry points are listed in [thermal project verification](THERMAL_PROJECT.md). Gate 2 needs clean installation and external evaluation. Gates 3 and 4 additionally require platform qualification, completed security and licensing reviews, acceptance testing and owner approval for distribution.
