# NOVA Development Roadmap

A milestone is complete only when its acceptance tests, verification evidence, and documentation are satisfied. A checked-in prototype is not automatically a completed milestone.

## M0 — Engineering foundation

**Implemented; CI verified on Windows and Linux.**

- C++20 and CMake project.
- CSR sparse matrix.
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
- Robin transfer integration on validated exterior edges.
- Mixed Dirichlet-Neumann and Dirichlet-Robin boundary benchmarks.
- Triangle orientation-independent integration.
- Harmonic patch and manufactured quartic-solution benchmarks.
- L2, H1-seminorm, and H1 error calculations using Duffy-transformed Gaussian quadrature.
- Nodal and integral-norm mesh-refinement regression tests.

Remaining:
- Reusable element and integration interfaces shared across physics.
- Explicit general H1 function-space abstraction and DOF manager.
- Pure-Neumann nullspace handling.
- P2 triangle elements.
- Mesh refinement and element-quality metrics.
- Generic residual and tangent assembly interfaces.
- Independent benchmark cross-check and numerical evidence report.

## M2 — Solver infrastructure

**In progress.**

Implemented:
- Conjugate Gradient for symmetric positive-definite systems.
- Optional Jacobi diagonal preconditioning for CG.
- BiCGSTAB for nonsymmetric systems.
- Optional Jacobi preconditioning for BiCGSTAB.
- Restarted GMRES for nonsymmetric systems.
- Right Jacobi preconditioning for GMRES.
- Modified Gram-Schmidt Arnoldi orthogonalization with a second reorthogonalization pass.
- Givens rotations for the GMRES least-squares problem.
- True-residual recomputation before convergence is accepted.
- Input validation, breakdown diagnostics, and iteration-limit tests.

Remaining:
- Incomplete Cholesky preconditioning.
- Sparse matrix block operations.
- Newton nonlinear solver.
- Line search and damping.
- Residual and Jacobian diagnostics.
- Solver configuration and reproducible logs.

## M3 — General scalar physics
- Thermal conduction.
- Volumetric and boundary sources.
- Temperature-dependent material properties.
- Stationary and transient scalar diffusion.
- Parameter sweeps.

## M4 — Structural mechanics
- Vector-valued H1 spaces.
- Linear elasticity.
- Plane stress and plane strain.
- 3D elasticity.
- Stress and strain recovery.
- Rigid-body and patch tests.

## M5 — Electromagnetics and optics
- Complex-valued algebra.
- H(curl) function spaces.
- Nedelec elements.
- Frequency-domain Maxwell equations.
- Eigenvalue formulations.
- Boundary conditions and absorbing layers.

## M6 — Multiphysics framework
- Shared model and field registry.
- Segregated coupling.
- Monolithic block assembly.
- Nonlinear material coupling.
- Coupled electrothermal and thermoelastic examples.

## M7 — Geometry, meshing, and visualization
- CAD import and geometry model.
- Mesh generation and quality inspection.
- 2D and 3D results visualization.
- Derived quantities and export.
- Project persistence.

## M8 — Desktop and commercial platform
- Qt desktop application.
- Model Builder.
- Materials and physics browsers.
- Study and solver configuration.
- Results viewer.
- Python automation.
- Licensing, installer, update, diagnostics, and support workflows.

## M9 — Scale and advanced capabilities
- Parallel assembly and solution.
- Adaptive mesh refinement.
- Domain decomposition.
- GPU acceleration where justified.
- Optimization and sensitivity analysis.
- Remote and batch execution.
