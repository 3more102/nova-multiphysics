# NOVA Architecture

## 1. Architectural objectives

NOVA is intended to become a commercial, general-purpose multiphysics CAE application. The architecture must support independent evolution of the desktop application, model representation, geometry and mesh services, numerical engine, physics modules, and results system.

The current implementation is an early numerical prototype. Its 1D and 2D scalar FEM cases are verification steps, not evidence of production readiness or broad multiphysics capability.

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

The current prototype includes a validated 2D triangle mesh, a structured rectangular-grid generator, and boundary-edge topology for structured rectangles. It does not yet include general CAD geometry, unstructured meshing, general physical groups, or adaptive mesh refinement.

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

Initial scalar modules use nodal P1 H1-like fields. The current 2D diffusion implementation is a focused triangle-P1 path with Dirichlet, Neumann, and Robin boundary contributions, not yet a general function-space framework. Electromagnetics will require compatible vector spaces such as Nedelec elements; mixed formulations must be supported explicitly rather than forced into a scalar interface.

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
- A 1D P1 finite-element Poisson benchmark.
- A validated 2D triangle mesh and structured rectangular mesh generator.
- P1 reference-triangle basis functions and positive degree-four quadrature.
- A 2D isotropic scalar diffusion assembler with spatially varying scalar diffusivity, volumetric sources, strong nodal Dirichlet constraints, and Neumann/Robin edge terms.
- Harmonic patch, mixed-boundary, and manufactured-solution benchmarks.
- L2 and H1 error-norm calculations for nodal P1 fields.
- Automated numerical and input-validation tests.

The 2D diffusion path currently lacks pure Neumann nullspace handling, anisotropic tensors, P2 elements, a general DOF manager, and adaptive meshing. These are explicit follow-up tasks.

## 4. Dependency policy
Every dependency must be evaluated for:
- Commercial redistribution rights.
- License compatibility and notice obligations.
- Security and maintenance history.
- Supported platforms and compiler versions.
- ABI stability and integration cost.

The core FEM algorithms must not silently depend on a third-party multiphysics solver.
