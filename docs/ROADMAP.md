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
- Per-component pure-Neumann compatibility and zero-area-mean gauge for P1/P2 triangles and P1 tetrahedra.
- Rejection of folded shared-edge neighbors and unused triangle-mesh vertices.
- Triangle orientation-independent integration.
- Orientation-aware 2D triangle quality analysis (minimum angle, mean ratio, element area and longest edge).
- 3D tetrahedral mesh quality analysis (total volume, mean ratio and minimum dihedral angle), exposed through pre-solve CLI inspection and user-selected quality gates.
- Conforming uniform 1-to-4 refinement with shared edge midpoints, exterior boundary splitting, parent-element maps and P1 nodal prolongation.
- Harmonic patch, P2 exact-quadratic, and manufactured quartic-solution benchmarks.
- L2, H1-seminorm, and H1 error calculations using Duffy-transformed Gaussian quadrature.
- Nonlinear 1D P1 diffusion residual and consistent Jacobian assembly.
- Manufactured nonlinear diffusion benchmark with Newton solution and mesh-refinement regression tests.

Remaining:
- Reusable element and integration interfaces shared across physics.
- Local/graded adaptive mesh refinement, error indicators, quality-based marking and solver mesh-transfer workflows.
- Generic residual and tangent assembly interfaces.
- Independent benchmark cross-check and numerical evidence report.

## M2 — Solver infrastructure

**In progress.**

Implemented:
- Conjugate Gradient for symmetric positive-definite systems.
- Optional Jacobi diagonal preconditioning for CG.
- Natural-order zero-fill incomplete Cholesky (IC(0)) preconditioning for CG, with diagonal equilibration and explicit rejection of unresolved pivots.
- BiCGSTAB for nonsymmetric systems.
- Optional Jacobi preconditioning for BiCGSTAB.
- Restarted GMRES with modified Gram-Schmidt reorthogonalization and Givens rotations.
- Optional right Jacobi preconditioning for GMRES.
- Newton nonlinear solver with GMRES or BiCGSTAB linear subsolves.
- Armijo backtracking line search, including rejection of non-finite trial residuals.
- True-residual recomputation before linear convergence is accepted.
- Input validation, breakdown diagnostics, and nonlinear regression tests.

Remaining:
- Sparse matrix block operations and multilevel preconditioners.
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
- 3D P1 tetrahedral adaptive Backward Euler/Crank–Nicolson with accepted-state streaming, once-per-study static FEM preparation, separate full/half system caches and exact-matrix IC(0) factor reuse across repeated time steps.
- 2D P1/P2 adaptive step-doubling prepares spatial operators once per study, keeps separate full/half dt-dependent base operators, and reuses IC(0) only for identical constrained CSR matrices while reevaluating time-dependent loads and boundary data at each time level.
- Exact H1 nodal DOF restarts with nonzero physical start times.
- Spatially varying positive volumetric heat capacity and conductivity; prescribed time-varying scalar/tensor conductivity and volumetric capacity in the 3D transient heat path.
- Time-dependent volumetric heating, essential Dirichlet temperatures, Neumann outward heat fluxes, and Robin ambient convection.
- Validated, mutually exclusive exterior-edge boundary assignment; other edges insulated by default.
- Integrated discrete heat budget with Dirichlet reaction contributions.
- Manufactured-solution temporal convergence and source/boundary energy-balance tests.
- Cooling plate, convective boundary-flux, and adaptive time-step example executables.

Remaining:
- Nonlinear temperature-dependent constitutive properties; time-dependent material support beyond the 3D transient heat path.
- Broader global temporal error estimation for adaptive step-doubling.
- Independently recovered normal heat fluxes and physical energy-balance cross-checks.
- Large-scale and distorted-mesh verification against independent solvers.
- Parameter sweeps and model-level thermal study configuration.


### 3D time-dependent conduction (experimental C++ API)

Implemented:
- P1 tetrahedral Backward Euler and Crank–Nicolson with a consistent thermal mass matrix.
- Optional additive time-dependent volumetric heating, with Backward Euler endpoint sampling and Crank–Nicolson trapezoidal load integration in fixed and adaptive studies.
- Fixed-mesh adaptive step doubling for both schemes, with full-step/two-half-step error control, accepted/rejected attempt history and an accepted-state observer.
- Versioned, mesh-validated text checkpoints for exact 3D transient nodal-state restart; physics model reconstruction remains caller-owned.
- Once-per-study static FEM preparation plus separate full-step and half-step dt-dependent operator caches; exact CSR IC(0) factor reuse is observable through result counters.
- Scalar and SPD tensor conductivity, per-cell material tensors, and spatial or cellwise heat capacity.
- Dirichlet, Neumann and Robin boundaries with per-step integrated energy diagnostics.
- Streaming state observers, optional bounded-memory result history, and validated ParaView PVD/VTK output.
- Optional spatial source plus additive cellwise source; absent source fields represent unforced conduction.
- Regression tests for heterogeneous cellwise SPD tensors, cell-only heating, zero-source equilibrium and energy-budget closure.
- Instance-owned IC(0) factor reuse across identical fixed-step 2D and 3D transient matrices; exact coefficient or sparsity changes trigger factor rebuilds.
- Analytical single-mode time-step refinement verifies first-order Backward Euler and second-order Crank–Nicolson with nonzero start time and per-step energy closure.
- Gmsh MSH 2.2 and 4.1 ASCII tagged transient CLI for first-order tetrahedra, physical-group material/source/boundary mapping, and streamed PVD/VTK output with independent Python verification.
- Time-dependent 3D Dirichlet nodal/face values, Neumann face fluxes, and Robin transfer/ambient data, with Backward Euler endpoint sampling and Crank–Nicolson endpoint operator/load integration.
- Time-dependent 3D scalar/tensor conductivity and volumetric capacity, with conservative enthalpy storage, endpoint-consistent Backward Euler/Crank–Nicolson operators, adaptive step-doubling support and energy-budget closure.

Remaining:
- Spatial error-driven refinement coupled to local tetrahedral refinement and mesh transfer; broader adaptive temporal-error estimators.
- Gmsh binary and higher-order element import, and independent cross-solver thermal benchmarks.
- Coupling to a reusable adaptive spatial-error estimator and mesh-transfer workflow.

### 3D steady conduction flux recovery (experimental)

Implemented:
- Element-averaged P1 conductive flux recovery for scalar and SPD tensor materials.
- Optional cellwise-conservative integrated heat rates on shared tetrahedral faces.
- Equal-and-opposite interior face rates, prescribed Neumann/Robin rates, and Dirichlet reaction outlets.
- Area-weighted minimum-change reconciliation with component-scaled sparse solves and local per-cell energy-budget checks.
- Raw P1 interior normal-flux continuity diagnostics on all shared tetrahedral faces, with per-face jumps, per-cell weighted indicators and a global indicator (not a full residual estimator).

Remaining:
- A complete residual-based 3D estimator including cell-source and exterior-boundary residuals, coupled to graded local tetrahedral refinement.
- H(div)-conforming vector flux spaces and enforced pointwise normal-flux continuity.
- Independent distorted-mesh and multi-material flux-convergence benchmarks.
- Export of recovered face rates to result files and visualization.

## M4 — Structural mechanics
- Vector-valued H1 spaces.
- Linear elasticity.
- Plane stress and plane strain.
- 3D elasticity.
- Stress and strain recovery.
- Rigid-body and patch tests.

## M5 — Electromagnetics and optics
Implemented (experimental):
- Scalar electrostatic Poisson equation on P1 tetrahedra with scalar or SPD tensor permittivity.
- Spatial and additive per-cell free charge, fixed potentials, prescribed normal displacement, pure-Neumann gauge and stored-field-energy diagnostics.
- Fixed-potential electrode charge diagnostics on validated exterior faces, using winding-independent normals and explicitly documented element-averaged displacement.
- Experimental multi-terminal Maxwell capacitance matrix extraction through unit-voltage basis solves, with an explicit reciprocity-defect diagnostic.

Remaining:
- Floating conductor potential unknowns from prescribed total charge and discrete electrode-reaction charge extraction.
- Locally equilibrated H(div) displacement recovery and Gauss-law error estimation.
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

Implemented (experimental):
- Gmsh MSH 2.2 and 4.1 ASCII first-order tetrahedral import with physical-group names, sparse IDs, exterior-face tagging and shared topology validation.
- Pre-solve tetrahedral mesh-quality reporting and user-selected minimum-mean-ratio screening gates for steady and transient Gmsh studies.

Remaining:
- CAD import and geometry model.
- Mesh generation and quality-driven remeshing.
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
