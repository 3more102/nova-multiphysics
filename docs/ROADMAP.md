# NOVA Development Roadmap

A milestone is complete only when its acceptance tests, verification evidence, and documentation are satisfied. A checked-in prototype is not automatically a completed milestone.

## M0 — Engineering foundation

**Implemented; CI verified on Windows and Linux.**

- C++20 and CMake project.
- CSR sparse matrix and O(nnz) scaled sparse linear combinations.
- 1D P1 Poisson FEM benchmark.
- Automated tests and continuous integration.

## M1 — FEM kernel

**In progress.**

Implemented:
- Triangle mesh data model with scale-aware geometric validation.
- Structured rectangular mesh generator with two P1 triangles per cell.
- Linear P1 and quadratic P2 reference-triangle H1 basis functions.
- Positive six-point degree-four triangle quadrature.
- Jacobian-based physical basis gradients.
- Consistent P1/P2 scalar H1 mass and stiffness matrix assembly for transient physics.
- Shared continuous H1 P1/P2 2D scalar isotropic diffusion assembly with orientation-independent physical gradients.
- Spatially varying positive scalar diffusivity and volumetric source.
- Strong symmetric Dirichlet elimination for vertex DOFs and P2 boundary-edge midpoint DOFs.
- Three-point Gauss Neumann outward-flux integration on validated exterior edges for P1/P2.
- Three-point Gauss Robin transfer integration on validated exterior edges for P1/P2.
- Mixed Dirichlet-Neumann and Dirichlet-Robin boundary benchmarks.
- Triangle orientation-independent integration.
- Orientation-aware 2D triangle quality analysis (minimum angle, mean ratio, element area and longest edge).
- Conforming uniform 1-to-4 refinement with shared edge midpoints, exterior boundary splitting, and parent-element maps.
- P1 and P2 solution prolongation across uniform refinement, with P2 transfer preserving the coarse quadratic field at all refined H1 DOFs.
- P1 residual-based cell and interior-flux-jump error indicators for positive elementwise-constant isotropic diffusivity, plus deterministic Dörfler bulk marking.
- Harmonic patch, P2 exact-quadratic, and manufactured quartic-solution benchmarks.
- L2, H1-seminorm, and H1 error calculations using Duffy-transformed Gaussian quadrature.
- Nonlinear 1D P1 diffusion residual and consistent Jacobian assembly.
- Manufactured nonlinear diffusion benchmark with Newton solution and mesh-refinement regression tests.

Remaining:
- Reusable element and integration interfaces shared across physics.
- Pure-Neumann nullspace handling.
- Conforming local/graded adaptive mesh refinement, physical-group transfer, and an automated solve-estimate-mark-refine workflow.
- Generic residual and tangent assembly interfaces.
- Independent benchmark cross-check and numerical evidence report.

## M2 — Solver infrastructure

**In progress.**

Implemented:
- Conjugate Gradient for symmetric positive-definite systems.
- Optional Jacobi diagonal preconditioning for CG.
- BiCGSTAB for nonsymmetric systems.
- Optional Jacobi preconditioning for BiCGSTAB.
- Restarted GMRES with modified Gram-Schmidt reorthogonalization and Givens rotations.
- Optional right Jacobi preconditioning for GMRES.
- Newton nonlinear solver with GMRES or BiCGSTAB linear subsolves.
- Armijo backtracking line search, including rejection of non-finite trial residuals.
- True-residual recomputation before linear convergence is accepted.
- Input validation, breakdown diagnostics, and nonlinear regression tests.

Remaining:
- Incomplete Cholesky preconditioning.
- Sparse matrix block operations.
- Inexact Newton forcing terms and nonlinear convergence diagnostics.
- Continuation and parameter-dependent nonlinear studies.
- Solver configuration and reproducible logs.

## M3 — General scalar physics

**In progress: transient thermal prototypes with P1/P2 elements implemented.**

Implemented:
- Stationary scalar conduction/diffusion with P1/P2 H1 elements.
- Consistent scalar capacity and stiffness FEM operators.
- Backward Euler and Crank–Nicolson implicit transient time integration.
- Adaptive, fixed-mesh time stepping by step doubling, absolute/relative local error control, accepted/rejected attempt histories, minimum/maximum step size and attempt limits.
- Exact H1 nodal DOF restarts with nonzero physical start times.
- Spatially varying positive (time-independent) volumetric heat capacity and conductivity.
- Time-dependent volumetric heating, essential Dirichlet temperatures, Neumann outward heat fluxes, and Robin ambient convection.
- Validated, mutually exclusive exterior-edge boundary assignment; other edges insulated by default.
- Integrated discrete heat budget with Dirichlet reaction contributions.
- Manufactured-solution temporal convergence and source/boundary energy-balance tests.
- Cooling plate, convective boundary-flux, and adaptive time-step example executables.

Remaining:
- Nonlinear temperature-dependent constitutive properties and time-dependent materials.
- Optimized adaptive stepping with shared preassembled FEM operators, factorization/preconditioner reuse and a global temporal error estimator.
- Independently recovered normal heat fluxes and physical energy-balance cross-checks.
- Large-scale and distorted-mesh verification against independent solvers.
- Parameter sweeps and model-level thermal study configuration.

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
