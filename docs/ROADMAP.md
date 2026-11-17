# NOVA Development Roadmap

A milestone is complete only when its acceptance tests, verification evidence, and documentation are satisfied. A checked-in prototype is not automatically a completed milestone.

## M0 — Engineering foundation

**Implemented; CI verified on Windows and Linux.**

- C++20 and CMake project.
- CSR sparse matrix.
- Conjugate-gradient solver.
- 1D P1 Poisson FEM benchmark.
- Automated tests and continuous integration.

## M1 — FEM kernel

**In progress.**

Implemented:
- Triangle mesh data model with scale-aware geometric validation.
- Structured rectangular mesh generator with two P1 triangles per cell.
- Linear P1 reference-triangle basis functions.
- Positive six-point degree-four triangle quadrature.
- Jacobian-based physical basis gradients.
- 2D scalar isotropic diffusion assembly.
- Spatially varying positive scalar diffusivity and volumetric source.
- Strong nodal Dirichlet elimination.
- Neumann outward-flux integration on validated exterior edges.
- Mixed Dirichlet-Neumann boundary benchmark.
- Triangle orientation-independent integration.
- Harmonic patch and manufactured quartic-solution benchmarks.
- L2, H1-seminorm, and H1 error calculations using Duffy-transformed Gaussian quadrature.
- Nodal and integral-norm mesh-refinement regression tests.

Remaining:
- Reusable element and integration interfaces shared across physics.
- Explicit general H1 function-space abstraction and DOF manager.
- Robin boundary conditions and pure-Neumann nullspace handling.
- P2 triangle elements.
- Mesh refinement and element-quality metrics.
- Generic residual and tangent assembly interfaces.
- Independent benchmark cross-check and numerical evidence report.

Acceptance: verified 1D and 2D scalar diffusion with mesh-convergence studies, boundary-condition coverage, and reusable function-space/assembly interfaces.

## M2 — Solver infrastructure

**In progress.**

Implemented:
- Unpreconditioned Conjugate Gradient.
- Optional Jacobi diagonal preconditioning.
- Explicit validation of Jacobi diagonal entries.
- True-residual recomputation before convergence is accepted.
- Regression tests on a scaled SPD diagonal system.

Remaining:
- Incomplete Cholesky preconditioning.
- GMRES and BiCGSTAB for nonsymmetric systems.
- Sparse matrix block operations.
- Newton nonlinear solver.
- Line search and damping.
- Residual and Jacobian diagnostics.
- Solver configuration and reproducible logs.

Acceptance: verified linear and nonlinear benchmark suite with documented convergence and breakdown behavior.

## M3 — General scalar physics
- Thermal conduction.
- Volumetric and boundary sources.
- Temperature-dependent material properties.
- Stationary and transient scalar diffusion.
- Parameter sweeps.

Acceptance: energy balance, manufactured solutions, and transient convergence tests.

## M4 — Structural mechanics
- Vector-valued H1 spaces.
- Linear elasticity.
- Plane stress and plane strain.
- 3D elasticity.
- Stress and strain recovery.
- Rigid-body and patch tests.

Acceptance: canonical elasticity benchmarks and reaction-force balance.

## M5 — Electromagnetics and optics
- Complex-valued algebra.
- H(curl) function spaces.
- Nedelec elements.
- Frequency-domain Maxwell equations.
- Eigenvalue formulations.
- Boundary conditions and absorbing layers.

Acceptance: waveguide and cavity benchmarks with mesh-convergence evidence.

## M6 — Multiphysics framework
- Shared model and field registry.
- Segregated coupling.
- Monolithic block assembly.
- Nonlinear material coupling.
- Coupled electrothermal and thermoelastic examples.

Acceptance: coupled benchmarks with independent energy and residual checks.

## M7 — Geometry, meshing, and visualization
- CAD import and geometry model.
- Mesh generation and quality inspection.
- 2D and 3D results visualization.
- Derived quantities and export.
- Project persistence.

Acceptance: end-to-end desktop-independent workflow for representative models.

## M8 — Desktop and commercial platform
- Qt desktop application.
- Model Builder.
- Materials and physics browsers.
- Study and solver configuration.
- Results viewer.
- Python automation.
- Licensing, installer, update, diagnostics, and support workflows.

Acceptance: usability-tested release candidate with documented validation scope.

## M9 — Scale and advanced capabilities
- Parallel assembly and solution.
- Adaptive mesh refinement.
- Domain decomposition.
- GPU acceleration where justified.
- Optimization and sensitivity analysis.
- Remote and batch execution.

Acceptance: published scaling, accuracy, and reproducibility reports.
