# Release acceptance criteria

Current status is an engineering prototype. The owner approves every release against an identified commit, reproducible build artifacts and recorded verification evidence.

| Target | Evidence required to close the gate | Open work |
| --- | --- | --- |
| 0.9 · engineering alpha | Successful Windows/Linux builds of the integrated revision. Reproducible 2D and 3D scalar and heat examples. Analytical patch, spatial/temporal convergence, orientation, boundary, conservation and failure-mode tests. Published numerical limits and performance baselines. | Integrate the dependency chain and qualify further use cases. The [3D thermal dossier](THERMAL_3D_VERIFICATION.md) records the first external FEM comparison, sanitizer/static-analysis evidence and single-host serial baselines. |
| 1.0 · thermal/scalar CAE | Stable public FEM and Python APIs, versioned project save/load and migration tests, a Qt workflow from geometry through results, reproducible installers, update/rollback, verified material and mesh import, and user/support documentation. Archive installation and customer-workflow tests for each supported OS. | Product workflow, installations, field validation, ownership review and licensing gates are open. |
| 2.0 · structural and coupled physics | Verified plane stress/strain and 3D elasticity, constitutive tangents, stress recovery, electrothermal and thermoelastic reference studies, adaptive transfer and conservation on shared field interfaces. | Structural and coupling modules and independent validation cases are open. |
| 3.0 · advanced CAE | Verified fluid, EM and acoustic studies, distributed reproducibility and scale measurements, with hardware/software-specific CPU and accelerator comparisons. | Advanced physics and industrial workloads are open. |

## Commercial and V&V evidence

- The owner records contributor provenance, rights to redistribute each shipped component and legal approval of the selected terms, EULA, notices and privacy policy. Each deliverable needs an SPDX-compatible SBOM, including runtime and installer dependencies.
- Product activation needs trials and tier/feature entitlements, renewal/recovery, offline and online flows, and tamper tests. Use a reviewed cryptographic library with public verification keys in the client and a separately controlled signing service. Enterprise floating licensing is a later gate.
- Each physics dossier identifies its intended use, mesh/time tolerances, analytical references, independent solver comparisons, experimental data and measurement uncertainty where available. Record quantitative acceptance bands and the limits of the validated domain.
- The V&V program references ASME V&V 20 for CFD/heat transfer and V&V 10 for computational solid mechanics. The evidence for this prototype is software verification. Experimental validation and standards assessments are open. [ASME VVUQ standards](https://www.asme.org/codes-standards/vvuq-standards)
