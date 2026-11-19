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

### Steady 3D conduction on a tagged solid

Import a conforming first-order tetrahedral Gmsh MSH 2.2 ASCII mesh, assign
conductivities by volume ID and fixed temperature, outward flux or convection
by exterior surface ID. The result file contains nodal temperature, volume IDs
and a heat-flux vector per tetrahedron. When region-specific heat sources are
provided, VTK also contains effective per-cell `heat_generation` (W/m³),
including any uniform `--source` contribution.

```sh
./build/nova_gmsh_diffusion_3d examples/gmsh_diffusion_3d/two_layer_wall.msh --material 11 2 --material 22 6 --temperature 101 400 --convection 102 12 300 --vtk wall.vtk
```

The included wall has a 340 K interface and transmits 240 W across its 1 m²
section. Assign distinct internal heat-generation densities by physical volume
with `--heat-source 11 120 --heat-source 22 360` (W/m³); signed
per-volume sources add to the optional uniform `--source Q`. The verification suite also measures temperature and flux convergence
under volumetric heating. [Formulation, input conventions and measurements](docs/STEADY_CONDUCTION_3D.md)

### Transient 3D heat on tetrahedral solids (experimental)

NOVA now solves `C(x) ∂T/∂t − div(K(x) grad T) = Q(x)` using a
**consistent P1 tetrahedral thermal mass matrix** and the existing 3D
conductivity integration kernel. Backward Euler and Crank–Nicolson time
stepping support scalar or SPD tensor conductivity, positive spatial or
per-cell heat capacity, volumetric heating, fixed temperature boundaries,
imposed outward flux and Robin convection. For now, coefficients and
boundary functions are **time-independent**.

The constant time-step CSR system is reused between steps; the solver
computes temperature corrections to avoid loss of precision for small
updates relative to absolute temperature. Each step reports heat content,
volumetric energy, imposed flux, convection, Dirichlet reactions and the
discrete energy-balance residual.

```sh
./build/nova_transient_heat_3d
```

See [3D transient formulation, API example and analytical tests](docs/TRANSIENT_HEAT_3D.md).
This is a C++ API and a structured-box demonstration, not yet an imported
Gmsh transient study with ParaView time-series output.

### Electrostatics 3D on tetrahedral solids (experimental)

The `nova/physics/electrostatics_3d.hpp` API reuses the independent P1 tetrahedral
Poisson engine to solve `-div(epsilon grad(V)) = rho_free`, with scalar or SPD
tensor permittivity, optional free-volume charge, prescribed electrode potentials
and exterior normal displacement. Unassigned faces have zero normal displacement.
It returns nodal potential, per-element electric field `E = -grad(V)`,
quadrature-averaged displacement `D = epsilon E`, and stored energy
`0.5 * integral(E dot D) dV` in SI units.

```sh
./build/nova_electrostatics_3d
```

The example is a 1 m cube with permittivity 2 F/m and two 10 V-separated
electrodes, with exact capacitor energy 100 J. Tests also cover anisotropic
tensor fields and rejection of invalid dielectric/charge inputs.
The API is experimental: no floating electrodes, nonlinear dielectric response
or time-dependent Maxwell equations. Per-cell displacement is not an
equilibrated face flux or local Gauss-law error estimator.

### Triangle mesh quality and uniform refinement

The mesh kernel now evaluates each triangle's signed area, minimum interior
angle, longest edge and mean-ratio shape quality (0 to 1; equilateral is 1).
It supports conforming uniform 1-to-4 refinement on a validated 2D mesh:
- Unique global edge midpoints shared between neighboring triangles.
- Four child triangles per parent, preserving orientation and domain area.
- Optional exterior-boundary-edge splitting that preserves edge direction.
- Parent-element indices and midpoint-parent edge records for mesh transfer.
- P1 nodal prolongation (coarse vertices copied; edge midpoint averages).
- Geometry and topology regression tests, including mixed element orientation.

API: `nova/fem/mesh_refinement_2d.hpp`. This is **uniform** refinement,
not error-driven adaptive mesh refinement. It does not yet preserve arbitrary
physical boundary groups or transfer P2 solution DOFs automatically.

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

Limitations: no temperature-dependent coefficients, adaptive mesh refinement, time integrators above order two, or independently reconstructed Dirichlet normal heat flux. Adaptive step doubling reassembles the FEM matrices for each trial and is not optimized for large systems. This prototype is not validated for engineering design use.

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

The 2D stationary solver requires an anchored problem and a scalar diffusivity. The 3D stationary solver also supports balanced floating components and spatial SPD conductivity tensors. Incomplete Cholesky, error-driven mesh adaptation and a general multiphysics field registry remain outside this mainline configuration. The nonlinear 1D example is a verification model, not yet a general nonlinear material framework.

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
./build/nova_transient_heat_3d
./build/nova_diffusion_2d
./build/nova_transient_heat_2d
./build/nova_transient_flux_2d
./build/nova_adaptive_heat_2d
./build/nova_nonlinear_diffusion_1d
~~~

On Windows, use the corresponding executable paths under the build directory (or the Release subdirectory for multi-configuration generators).

## Experimental C++ SDK package

NOVA's built core can be installed as a relocatable CMake package, consumed by a separate C++20 project with `find_package(NOVA 0.8 CONFIG REQUIRED)` and `target_link_libraries(... NOVA::Core)`.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cmake --install build --config Release --prefix ./nova-sdk
cmake -S tests/package_consumer -B build-package-consumer -DCMAKE_PREFIX_PATH="$(pwd)/nova-sdk"
cmake --build build-package-consumer --config Release --parallel
ctest --test-dir build-package-consumer --output-on-failure
```

See [installation and external consumer verification](docs/SDK_INSTALL.md) for cross-platform instructions and limitations. Installation support **does not** grant commercial redistribution rights or imply certification.

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
