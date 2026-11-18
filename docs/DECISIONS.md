# Engineering decisions

| Date | Decision | Reason and verification consequence |
| --- | --- | --- |
| 2026-11-18 | Continue the 3D path from the exact head of #16 on a dedicated branch. | Its complete upstream ancestry includes main. Isolate sibling features until their own contracts, merge plans and tests agree. Record the base SHAs in [project status](PROJECT_STATUS.md). |
| 2026-11-18 | Use dimensionless edge vectors and unit face normals for tetrahedral shape metrics. | Uniform scaling must preserve shape. A shared-face pair also requires opposite outward windings. The mesh validator owns these checks and rejects vertices without element support. |
| 2026-11-18 | Return the normalized physical field and preserve the original CG iterate as solver diagnostics. | The temporary component pins define the matrix whose residual CG reports. A one-tetrahedron analytical matrix verifies the reported residual independently. |
| 2026-11-18 | Prioritize fixed-step 3D thermal using the existing 2D theta time integrator. | This reuses the verified stationary spatial discretization and permits temporal accuracy, flux and conservation checks before adaptivity or nonlinear materials are added. |
| 2026-11-18 | Defer the distribution model and commercial terms to the owner and legal review. | Contributor rights, third-party notices and exact shipped binaries must be documented before a commercial release decision. |
| 2026-11-18 | Share the theta time integrator, energy budget and positive simplex quadrature across scalar studies. | This keeps the 2D/3D endpoint semantics identical and confines geometry to element and boundary assembly. Only time levels with nonzero weight evaluate load callbacks. |
| 2026-11-18 | Initialize constrained identity rows to the exact prescribed trace in CG. | A split run exposed a discrepancy above 1e-12 when the continuous run accepted a slightly approximate Dirichlet value and the restart projected it. Keeping identity-row updates zero makes uninterrupted and restarted solves agree. |
| 2026-11-18 | Keep the independent reference backend in an optional Python verification script. | scikit-fem assembles separate operators and SciPy uses sparse LU. The C++ deliverable has no dependency on either package. |
| 2026-11-18 | Use the same compatibility, gauge, normalization and symmetric elimination in 2D and 3D stationary scalar solves. | The element supplies its own basis integrals. P2 triangle vertex basis functions integrate to zero. The raw CG diagnostic and mean-adjusted physical field retain their separate meanings. |
