# NOVA Multiphysics

**Engineering the Physics of Tomorrow.**

NOVA is an independent C++20 multiphysics simulation platform under development. Its long-term objective is a commercial, general-purpose CAE environment for geometry, meshing, finite-element modeling, numerical solution, multiphysics coupling, and scientific visualization.

> **Development status:** Early engineering prototype. It is not production-certified and must not be used as the sole basis for engineering design decisions.

## Current numerical prototype

### 1D Poisson
- Uniform 1D mesh.
- Linear P1 finite elements.
- Two-point Gauss quadrature.
- Strong endpoint Dirichlet conditions.
- CSR assembly and conjugate-gradient solution.
- Analytical quadratic benchmark.

### 2D scalar diffusion

The prototype solves

\[
-\nabla\cdot(k\nabla u)=f \quad \text{in }\Omega
\]

on triangular P1 meshes using:
- Structured rectangular mesh generation, split into two triangles per cell.
- Positive six-point degree-four reference-triangle quadrature.
- Jacobian-based physical basis gradients.
- Positive scalar isotropic diffusivity, including spatially varying coefficients.
- Volumetric source integration.
- Strong nodal Dirichlet elimination.
- Neumann outward-flux integration on validated exterior edges.
- Mixed Dirichlet-Neumann boundary benchmark.
- Orientation-independent triangle handling and scale-aware degenerate-element rejection.

### Quantitative verification
- Linear harmonic patch test with nonzero boundary values.
- Variable diffusivity and source verification.
- Reversed triangle orientation and invalid-geometry tests.
- Manufactured quartic solution on 4×4, 8×8, and 16×16 meshes.
- Nodal maximum error, L2 error, H1 seminorm, and H1 norm.
- Duffy-transformed 5×5 Gauss-Legendre integration for error norms.

### Linear solver
- CSR sparse matrix-vector products.
- Conjugate Gradient for symmetric positive-definite systems.
- Optional Jacobi diagonal preconditioning.
- Validation of positive finite Jacobi diagonal entries.
- True-residual recomputation before reporting convergence.
- Explicit breakdown and convergence diagnostics.

Jacobi preconditioning is diagonal scaling. It is not a substitute for incomplete-factorization or multilevel preconditioners on difficult large-scale systems.

The 2D implementation does not yet support Robin conditions, pure Neumann nullspace handling, anisotropic diffusion tensors, P2 elements, adaptive meshing, or a general DOF manager.

## Build requirements
- CMake 3.20 or newer for command-line builds.
- CMake 3.21 or newer for the supplied presets.
- A C++20 compiler (GCC 11+, Clang 14+, or MSVC 2022 recommended).
- Ninja for the supplied presets.

### Preset workflow
~~~bash
cmake --preset default
cmake --build --preset default
ctest --preset default
~~~

### Generator-independent workflow
~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
~~~

Run the examples:
~~~bash
./build/nova_poisson_1d
./build/nova_diffusion_2d
~~~

On Windows, use the corresponding executable paths under the build directory (or the Release subdirectory for multi-configuration generators).

## Planned product architecture
- **Application:** Qt desktop interface and project/model management.
- **Model layer:** geometry entities, materials, fields, parameters, units, boundary conditions, and studies.
- **Geometry and mesh:** CAD exchange, topology, mesh generation, validation, and refinement.
- **Numerical core:** function spaces, element definitions, quadrature, DOF mapping, assembly, sparse algebra, and solvers.
- **Physics:** thermal, electrostatic, structural, electromagnetic, fluid, acoustic, and semiconductor interfaces.
- **Multiphysics:** segregated and monolithic coupling.
- **Results:** datasets, derived quantities, visualization, and reporting.
- **Automation:** Python API and reproducible job execution.

## Engineering principles
- Numerical correctness before feature breadth.
- Analytical and manufactured-solution verification.
- Explicit units, sign conventions, and boundary-condition semantics.
- Reproducible automated tests on Windows and Linux.
- No silent fallback when a solver fails.
- Third-party dependencies require commercial license and security review.

## Licensing
The project license has not yet been finalized. Do not redistribute this prototype commercially until ownership, dependency licenses, and distribution terms have been reviewed.

## Disclaimer
NOVA is experimental research and software engineering work. It is not a certified engineering tool. Users must independently verify all results before relying on them.
