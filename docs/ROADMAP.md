# NOVA Development Roadmap

This roadmap is capability-based. A milestone is complete only when its acceptance tests and documentation are satisfied.

## M0 — Engineering foundation

- C++20 and CMake project.
- CSR sparse matrix.
- Conjugate-gradient solver.
- 1D P1 Poisson FEM benchmark.
- Unit tests and continuous integration.

Acceptance: reproducible build and passing analytical benchmark tests on Windows and Linux.

## M1 — FEM kernel

- Mesh and element abstractions.
- Reference segment and triangle elements.
- Basis functions and derivatives.
- Gaussian quadrature.
- Jacobian mappings.
- H1 function spaces.
- DOF manager.
- Generic element and global assembly.
- Strong Dirichlet constraints.

Acceptance: verified 1D and 2D scalar diffusion with convergence studies.

## M2 — Solver infrastructure

- Preconditioned Krylov methods.
- Sparse matrix block operations.
- Newton nonlinear solver.
- Line search and damping.
- Residual and Jacobian diagnostics.
- Solver configuration and reproducible logs.

Acceptance: verified linear and nonlinear benchmark suite.

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
