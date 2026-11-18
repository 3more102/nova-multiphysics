# Engineering decisions

| Date | Decision | Reason and verification consequence |
| --- | --- | --- |
| 2026-11-18 | Continue the 3D path from the exact head of #16 on a dedicated branch. | Its complete upstream ancestry includes main. Isolate sibling features until their own contracts, merge plans and tests agree. Record the base SHAs in [project status](PROJECT_STATUS.md). |
| 2026-11-18 | Use dimensionless edge vectors and unit face normals for tetrahedral shape metrics. | Uniform scaling must preserve shape. A shared-face pair also requires opposite outward windings. The mesh validator owns these checks and rejects vertices without element support. |
| 2026-11-18 | Return the normalized physical field and preserve the original CG iterate as solver diagnostics. | The temporary component pins define the matrix whose residual CG reports. A one-tetrahedron analytical matrix verifies the reported residual independently. |
| 2026-11-18 | Prioritize fixed-step 3D thermal using the existing 2D theta time integrator. | This reuses the verified stationary spatial discretization and permits temporal accuracy, flux and conservation checks before adaptivity or nonlinear materials are added. |
| 2026-11-18 | Defer the distribution model and commercial terms to the owner and legal review. | Contributor rights, third-party notices and exact shipped binaries must be documented before a commercial release decision. |
