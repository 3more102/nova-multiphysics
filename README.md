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

on continuous triangular P1 and P2 meshes using:
- Structured rectangular mesh generation, split into two triangles per cell.
- Positive six-point degree-four reference-triangle quadrature.
- Jacobian-based physical basis gradients.
- Positive scalar isotropic diffusivity, including spatially varying coefficients.
- Volumetric source integration.
- Strong symmetric Dirichlet elimination on vertex and (for P2) edge-midpoint degrees of freedom.
- Neumann outward-flux integration on validated exterior edges.
- Robin boundary integration with non-negative transfer coefficients and ambient values.
- Mixed Dirichlet-Neumann and Dirichlet-Robin benchmarks.
- L2 and H1 error norms with Duffy-transformed Gaussian quadrature.
- Orientation-independent triangle handling and scale-aware degenerate-element rejection.

### Transient 2D heat diffusion
The continuous P1/P2 H1 kernel solves

\[
\rho c\frac{\partial T}{\partial t}-\nabla\cdot(k\nabla T)=Q(x,y,t)
\]

using:
- Consistent mass/capacity and stiffness matrices in sparse CSR format.
- Positive spatially varying, time-independent volumetric heat capacity and conductivity.
- Backward Euler (first-order) or Crank–Nicolson (second-order for sufficiently smooth solutions) time integration.
- Adaptive time stepping with two-half-step error estimates, automatic rejection, absolute/relative RMS DOF tolerances and configurable step limits.
- Restart from an exact P1/P2 nodal DOF vector at a user-defined physical start time.
- Time-dependent volume sources and prescribed Dirichlet temperature traces, including P2 edge midpoint DOFs.
- Time-dependent Neumann outward heat flux and Robin convection with ambient temperature.
- Homogeneous natural (zero-flux) conditions on unassigned edges.
- Integrated energy diagnostics: heat content, volumetric energy, Neumann/Robin exchange, and constrained Dirichlet reaction.
- Conjugate Gradient and per-step residual/convergence reports.
- Manufactured-solution convergence and energy-budget regression tests.

Examples: examples/transient_heat_2d/main.cpp (P2 cooling plate), examples/transient_flux_2d/main.cpp (P2 Crank–Nicolson boundary exchange), and examples/adaptive_heat_2d/main.cpp (error-controlled P2 Crank–Nicolson heating).

Limitations: the linear/adaptive thermal path has no temperature-dependent coefficients, adaptive mesh refinement, time integrators above order two, or independently reconstructed Dirichlet normal heat flux. Adaptive step doubling reassembles the FEM matrices for each trial and is not optimized for large systems. This prototype is not validated for engineering design use.

### Temperature-dependent 2D heat conduction (experimental)
A separate P1/P2 Backward Euler solver freezes positive nonlinear conductivity
`k(x,y,T)` at each Picard iteration, then reuses the linear thermal assembler
and solver. It supports the same time-dependent sources and Dirichlet, Neumann,
and Robin boundary-condition inputs. Scaled DOF updates, under-relaxation and
an explicit iteration limit control convergence. It rejects invalid
conductivity and nonconvergence rather than returning partial results.

See `docs/NONLINEAR_THERMAL_PICARD.md` and
`examples/nonlinear_heat_2d/main.cpp`. This is not yet an assembly-native
Newton solver. Material evaluation currently performs a linear triangle
search at each quadrature request and is unsuitable for large meshes.

### Nonlinear 1D diffusion
The prototype also includes a nonlinear P1 finite-element benchmark:

\[
-\frac{d}{dx}\left[(1+u^2)\frac{du}{dx}\right]=f(x).
\]

- Nonlinear residual and consistent tangent/Jacobian assembly.
- Three-point Gauss quadrature.
- Strong endpoint Dirichlet constraints.
- Newton iterations with Armijo backtracking.
- Manufactured solution \(u(x)=\sin(\pi x)\).
- Mesh-refinement checks at 8, 16, and 32 elements.

### Linear and nonlinear solver capabilities
- CSR sparse matrix-vector products.
- Conjugate Gradient for symmetric positive-definite systems.
- BiCGSTAB and restarted GMRES for general nonsymmetric systems.
- Optional Jacobi diagonal preconditioning for all three iterative linear solvers.
- Newton nonlinear solver with GMRES or BiCGSTAB linear subsolves.
- Armijo backtracking, including rejection of non-finite trial residuals.
- Explicit breakdown and convergence reporting.
- True-residual recomputation before linear convergence is accepted.

Jacobi preconditioning is diagonal scaling. It is not a substitute for incomplete-factorization or multilevel preconditioners on difficult large-scale systems.

The current prototype does not yet support incomplete Cholesky, stationary pure-Neumann nullspace handling, anisotropic diffusion tensors, adaptive meshing, or a general multiphysics field registry. The nonlinear 1D example is a verification model, not yet a general nonlinear material framework.

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
./build/nova_transient_heat_2d
./build/nova_transient_flux_2d
./build/nova_adaptive_heat_2d
./build/nova_nonlinear_heat_2d
./build/nova_nonlinear_diffusion_1d
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
