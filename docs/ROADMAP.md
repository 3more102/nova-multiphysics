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
- Tetrahedral 3D mesh data model with scale-aware volume validation, manifold face checks, structured box generation and outward boundary-face extraction.
- Structured rectangular mesh generator with two P1 triangles per cell.
- Linear P1 and quadratic P2 reference-triangle H1 basis functions.
- 3D reference-tetrahedron P1 H1 basis, physical gradient transformation and continuous vertex DOF space.
- Positive four-point degree-two tetrahedron quadrature and three-point triangular-face quadrature.
- 3D P1 scalar diffusion assembly and CG solve with spatially varying positive isotropic diffusivity, volume sources, strong Dirichlet faces/vertices, Neumann fluxes, Robin face exchange, and component-wise pure-Neumann compatibility and zero-mean handling.
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
- 3D tetrahedron quality diagnostics (signed/absolute volume, mean ratio and minimum dihedral angle).
- Conforming uniform 1-to-4 refinement with shared edge midpoints, exterior boundary splitting, parent-element maps and P1 nodal prolongation.
- One-level conforming red-green adaptive refinement with unique edge midpoints, transition elements, parent-cell maps and directed boundary-edge remapping.
- Area-weighted H1 gradient recovery indicators and deterministic Dörfler bulk marking for P1/P2 fields.
- Exact P1/P2 H1 solution transfer to uniform and locally refined triangular meshes, including continuity checks across coarse-element interfaces.
- Uniform and selective transient heat remeshing with directed Dirichlet, Neumann and Robin edge splitting and optional P1/P2 restart-state transfer.
- Harmonic patch, P2 exact-quadratic, and manufactured quartic-solution benchmarks.
- End-to-end adaptive thermal example: estimate, Dörfler mark, red-green remesh, P2 restart transfer and re-solve.
- L2, H1-seminorm, and H1 error calculations using Duffy-transformed Gaussian quadrature.
- Nonlinear 1D P1 diffusion residual and consistent Jacobian assembly.
- Manufactured nonlinear diffusion benchmark with Newton solution and mesh-refinement regression tests.

Remaining:
- Reusable element and integration interfaces shared across physics.
- 2D pure-Neumann nullspace handling.
- Multi-level graded refinement closure, residual-based error indicators, quality-aware marking, arbitrary physical-group transfer and automated solve-estimate-mark-refine cycles.
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
- Stationary 2D scalar conduction/diffusion with P1/P2 H1 elements and 3D tetrahedral scalar diffusion with P1 H1 elements.
- Consistent scalar capacity and stiffness FEM operators, including 3D P1 tetrahedra.
- Fixed-step 3D transient heat with time-dependent Dirichlet, Neumann and Robin faces, restart states and a discrete heat budget. Analytical and independent FEM comparisons are recorded in [thermal verification](THERMAL_3D_VERIFICATION.md).
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
- Adaptive time stepping for 3D heat, nonlinear thermal properties and higher-order tetrahedral spaces.
- General material models and time-dependent volume coefficients.
- Optimized adaptive stepping with shared preassembled FEM operators, factorization/preconditioner reuse and a global temporal error estimator.
- General-purpose face flux reconstruction and broader physical validation. An analytical mixed-boundary heat budget is covered by the 3D regression suite.
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
