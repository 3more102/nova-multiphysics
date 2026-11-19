# NOVA Architecture

## 1. Architectural objectives

NOVA is intended to become a commercial, general-purpose multiphysics CAE application. The architecture must support independent evolution of the desktop application, model representation, geometry and mesh services, numerical engine, physics modules, and results system.

The current implementation is an early numerical prototype. The tested code covers scalar finite elements in one, two and three dimensions. See [project status](PROJECT_STATUS.md) for the revision and numerical evidence.

## 2. Target subsystem boundaries

### Desktop application

Responsibilities:
- Model Builder and project navigation.
- Geometry and mesh interaction.
- Material and physics configuration.
- Study and solver configuration.
- Results visualization and report generation.

The desktop application must not own numerical algorithms. It submits immutable or versioned model snapshots to a solver process.

### Model and project layer

Responsibilities:
- Persistent identifiers for model entities.
- Parameters, units, variables, and expressions.
- Geometry topology references.
- Material assignments.
- Physics interfaces and boundary conditions.
- Study definitions and solver settings.
- Project serialization and schema migration.

A geometry edit must not silently invalidate unrelated assignments. Entity references need persistent IDs and explicit remapping rules.

### Geometry and mesh services

Responsibilities:
- CAD import and construction.
- Topology and entity classification.
- Mesh generation and validation.
- Boundary and domain markers.
- Mesh refinement and adaptation.
- Mapping between geometric entities and mesh entities.

External CAD and mesh libraries may be used after commercial license review. NOVA must own its internal mesh representation and the contracts that connect mesh data to the numerical core.

The current prototype includes a validated 2D triangle mesh, a structured rectangular-grid generator, boundary-edge topology for structured rectangles, and a reusable continuous H1 function-space/DOF layer supporting P1 and P2 triangular Lagrange fields. The mesh layer also owns tetrahedral geometry, boundary face extraction, triangle red-green refinement and P1/P2 solution transfer. General CAD import, general physical groups and 3D adaptivity remain open.

### Numerical core

Responsibilities:
- Reference elements and basis functions.
- Function spaces and field layouts.
- Degrees-of-freedom numbering.
- Quadrature and coordinate mappings.
- Element residual and tangent evaluation.
- Global sparse assembly.
- Constraints and boundary-condition enforcement.
- Linear, nonlinear, transient, and eigenvalue solution workflows.
- Error estimation and convergence diagnostics.

The numerical core is a proprietary engineering focus. Its algorithms require documented assumptions and independent verification.

### Physics framework

Each physics interface declares its fields, variables, constitutive data, constraints, and contributions to residual and Jacobian assembly. Physics modules must not bypass shared DOF, units, material, or solver services.

Initial scalar modules use nodal H1 Lagrange fields. The 2D scalar diffusion implementation now assembles P1 and P2 elements through the shared H1 function space and DOF manager. Both orders support strong Dirichlet elimination and weak Neumann/Robin boundary contributions; P2 boundary traces constrain edge midpoint DOFs. The nonlinear 1D diffusion example demonstrates a temperature/field-dependent scalar constitutive coefficient through caller-defined residual and tangent assembly. It is a verification case, not a general material framework.

Electromagnetics will require compatible vector spaces such as Nedelec elements; mixed formulations must be supported explicitly rather than forced into a scalar interface.

### Multiphysics framework

Two coupling strategies are required:
- Segregated: solve one field or physics block at a time, with controlled field updates.
- Monolithic: assemble a global coupled residual and Jacobian and solve for all dependent fields together.

The coupling layer owns dependency graphs, block structure, convergence criteria, relaxation, and field transfer.

### Results and visualization

The results subsystem stores mesh identity, field values, units, study metadata, solver diagnostics, and derived quantities. Visualization must consume results without mutating solver state.

The current numerical prototype also provides elementwise L2 and H1 error measurements against analytical fields. These are verification utilities, not yet a general adaptive error-estimation framework.

### Solver process

The long-term desktop product should execute simulations in a separate process. A versioned job manifest and result format will allow cancellation, progress reporting, remote execution, and future HPC integration.

## 3. Current code boundaries

The prototype currently contains:
- CSR sparse matrix storage and matrix-vector multiplication.
- Conjugate Gradient with optional Jacobi diagonal preconditioning and true-residual convergence verification.
- BiCGSTAB and restarted GMRES for general nonsymmetric systems, with optional Jacobi diagonal preconditioning.
- A Newton nonlinear solver with GMRES or BiCGSTAB linear subsolves and Armijo backtracking.
- A 1D P1 finite-element Poisson benchmark.
- A validated 2D triangle mesh and structured rectangular mesh generator.
- Reusable H1 P1/P2 reference-triangle basis functions and positive degree-four quadrature.
- A shared H1 P1/P2 isotropic scalar diffusion assembler with spatially varying scalar diffusivity, volumetric sources, strong nodal and boundary-edge Dirichlet constraints, and P1/P2 Neumann/Robin edge terms.
- Reusable consistent capacity/mass and stiffness operators for P1/P2 triangles and P1 tetrahedra. The 3D stationary and transient assemblies share barycentric quadrature, physical point mapping, scalar/tensor stiffness integration and exterior-face assignment checks.
- A geometry-independent theta method in `src/fem/heat_time_integration.cpp`. The 2D P1/P2 and 3D P1 adapters supply material operators, time-level volume and face/edge loads, and prescribed traces. Both receive the same solver and energy reports.
- Fixed-mesh adaptive time-step controllers with full-step/two-half-step error estimation, accepted/rejected trial recording, and exact DOF restarts. The 2D reference path still reassembles fixed operators per trial. The 3D P1 path prepares the material mass/stiffness operators and element Jacobians once, then reuses them across full, half and rejected trials.
- A 1D nonlinear P1 diffusion benchmark with a consistent tangent/Jacobian for (k(u)=1+u^2).
- Harmonic patch, mixed-boundary, manufactured-solution, and nonlinear diffusion benchmarks.
- L2 and H1 error-norm calculations for nodal P1 fields.
- Automated numerical and input-validation tests.

The 2D and 3D stationary diffusion adapters share a component nullspace implementation. The DOF map defines connectivity. Each adapter supplies its exact basis integrals and quadrature load magnitudes. The shared symmetric elimination compacts free triplets in place, and CG starts all identity rows at their prescribed values. The tetrahedral adapters accept scalar or symmetric positive-definite tensor conductivity. A generalized physics residual/tangent interface remains open. The nonlinear solver is not yet connected to a general material or multiphysics model.

## 4. Dependency policy

Every dependency must be evaluated for:
- Commercial redistribution rights.
- License compatibility and notice obligations.
- Security and maintenance history.
- Supported platforms and compiler versions.
- ABI stability and integration cost.

The core FEM algorithms must not silently depend on a third-party multiphysics solver.

## Module ownership and interfaces

| Layer | Owner and testable boundary |
| --- | --- |
| Mathematical core and solvers | `nova::linalg` owns CSR values and solve reports. `nova::nonlinear` consumes residual/Jacobian callbacks. Matrix and residual tests exercise these contracts. |
| Geometry and mesh | Meshes own coordinates/connectivity. The 3D validator verifies shared-face topology, orientation and supported scale. Future CAD entities need versioned IDs and checked remapping. |
| FEM | Function spaces borrow a mesh for their lifetime and own DOF maps. Assembly produces owned operators, loads and physical fields. Stationary 3D field values and gauged CG diagnostics have separate consumers. |
| Materials and physics | Callers currently supply coefficient and boundary callbacks. Future material records need units, ranges and provenance. Modules will contribute to shared field equations. |
| Coupling and studies | Planned studies own snapshots, field layouts, time settings and solver parameters. Coupling owns transfer, iteration state and convergence. Replay must reconstruct the same inputs. |
| Results and visualization | Results own field values and diagnostic metadata. Viewers and exporters consume immutable results tied to mesh and time. |
| Desktop and Python | Planned Qt and Python clients submit the same versioned study description to the engine. Verify round trips, cancellation and stable public entry points. |
| Commercial infrastructure | An entitlement service owns product feature grants. Its public verification interface belongs at application entry points. Numerical assembly has no dependency on license state. |


## Fixed-step thermal data flow

A fixed-step solver call borrows the problem and owns its assembled CSR matrices. The 3D adapter stores one absolute element Jacobian for reuse in time-dependent volume loads. The adaptive 3D controller extends this lifetime across all full-step, half-step and rejected trials: static mass/stiffness operators and Jacobians are assembled once, while time-step-dependent system matrices and time-dependent source/face loads are rebuilt only when required. The time integrator retains the current mass action and previous heat rate only when the scheme needs it. The matrix with fixed Dirichlet topology is reused on steps with zero Robin matrix. The mesh, callbacks and initial DOF state are never modified.

Both dimensional APIs expose their existing dimensional names for the common operator, diagnostic and result types. Results retain one CG solution vector per accepted step, so solution-history memory grows with DOFs times steps. Streaming result storage is an open design requirement.
