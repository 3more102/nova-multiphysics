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
- Per-connected-component pure-Neumann compatibility checks and zero-area-mean P1/P2 gauges, including mixed floating, Robin-anchored and Dirichlet-anchored components.
- L2 and H1 error norms with Duffy-transformed Gaussian quadrature.
- Orientation-independent triangle handling and scale-aware degenerate-element rejection.

### Steady 3D conduction on a tagged solid

Import a conforming first-order tetrahedral Gmsh MSH 2.2 or 4.1 ASCII mesh, assign
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

Use `--jump-report` on a solved Gmsh case to inspect raw P1 normal heat-flux
jumps across shared tetrahedral faces. This is a flux-continuity diagnostic,
not a complete PDE error estimator or a local-conservation certificate.

Inspect a mesh before simulation, or apply a user-selected quality gate before
the FEM solve. NOVA reports total tetrahedral volume, minimum mean ratio,
minimum dihedral angle and the worst element index:

```sh
./build/nova_gmsh_diffusion_3d tests/data/unit_tetrahedron.msh --check-mesh --min-mean-ratio 0.8
```

Add `--quality-report --min-mean-ratio q` to a thermal solve to reject a mesh
before assembly when its minimum mean ratio is below `q`. This is a geometric
screening threshold, not proof of solution accuracy or engineering certification.

### Transient 3D heat on tetrahedral solids (experimental)

NOVA now solves `C(x) ∂T/∂t − div(K(x) grad T) = Q_static(x) + Q_time(x,t)` using a
**consistent P1 tetrahedral thermal mass matrix** and the existing 3D
conductivity integration kernel. Backward Euler and Crank–Nicolson time
stepping support scalar or SPD tensor conductivity, positive spatial or
per-cell heat capacity, volumetric heating, fixed temperature boundaries,
imposed outward flux and Robin convection. An optional additive time-dependent volumetric source is sampled at the new
step time for Backward Euler and averaged at both endpoints for Crank–Nicolson.
Time-dependent Dirichlet nodal/face values, Neumann face fluxes and Robin
transfer/ambient data are also supported. Dirichlet values are imposed at each
endpoint; Neumann fluxes and Robin data use endpoint sampling for Backward
Euler and trapezoidal integration for Crank–Nicolson. A changing Robin
transfer coefficient changes the time-level operator and triggers exact CSR
cache invalidation; ambient-only changes affect the load. Conductivity and
capacity remain **time-independent**.

The constrained CSR system is reused while the time step and Robin
transfer matrix are unchanged; a changed transfer matrix rebuilds the affected
operator and, when necessary, its IC(0) factor. The solver computes temperature
corrections to avoid loss of precision for small updates relative to absolute
temperature. Each step reports heat content,
volumetric energy, imposed flux, convection, Dirichlet reactions and the
discrete energy-balance residual.

```sh
./build/nova_transient_heat_3d
./build/nova_transient_heat_3d --vtk-series thermal_case.pvd
```

The second command runs the same single 3D FEM solve and streams its
Dirichlet-projected initial field plus every completed step to separate
legacy tetrahedral `.vtk` frames. Open `thermal_case.pvd` in ParaView to
animate temperature over physical time. A validated PVD writer and an optional
per-step solver observer avoid repeated assembly/solves and avoid retaining
all intermediate nodal fields when `retain_step_solutions=false`.
See [3D transient formulation, API example and analytical tests](docs/TRANSIENT_HEAT_3D.md).

Experimental C++ users can also serialize and restore the exact nodal state at a
nonzero physical time with a [versioned, mesh-validated 3D thermal checkpoint](docs/THERMAL_CHECKPOINT_3D.md).
This is a state-only checkpoint; the same material model, boundary conditions
and numerical settings must be reconstructed before continuation.

Experimental **adaptive 3D timestep control** now compares one full step
with two half steps for Backward Euler or Crank–Nicolson. The accepted
half-step solutions preserve energy diagnostics and support restart. Its
`on_step` observer emits the Dirichlet-projected initial field and accepted
macro-step fields only; rejected and coarse full-step trial fields are not
published.
The adaptive study prepares static tetrahedral material and load operators once,
then keeps separate full-step and half-step system caches. Time-dependent Robin
transfer matrices are evaluated at each trial's physical endpoints; ambient
loads remain right-hand-side data. Exact constrained-CSR comparison reuses
operators and IC(0) factors when the matrix is unchanged, and rebuilds only
when the step size or Robin transfer matrix changes. The result exposes operator
and factor-build counters. See [adaptive 3D API and limitations](docs/ADAPTIVE_HEAT_3D.md).
The structured-box example is complemented by an imported-mesh CLI. It
reads Gmsh MSH 2.2 or 4.1 ASCII tetrahedra, maps conductivity, heat capacity, heat sources
and thermal boundary conditions from physical groups, and streams all solved
states to ParaView `.pvd` time series:

```sh
./build/nova_gmsh_transient_heat_3d examples/gmsh_diffusion_3d/two_layer_wall.msh \
  --material 11 2 --material 22 6 \
  --capacity 11 2000000 --capacity 22 3500000 \
  --initial 300 --dt 0.5 --steps 40 --scheme cn \
  --temperature 101 310 --convection 102 12 290 \
  --heat-source 22 120 --vtk-series heated_wall.pvd
```

The collection includes the projected initial field and every completed
time step, with nodal temperature and per-cell material, heat-generation and
heat-flux fields. See [Gmsh transient CLI and verification](docs/GMSH_TRANSIENT_HEAT_3D.md).

### Experimental 3D electrostatics on tetrahedral solids

The experimental P1 tetrahedral electrostatics API solves

\[
-\nabla\cdot(\epsilon\nabla V)=\rho_{free},
\qquad E=-\nabla V,
\qquad D=\epsilon E.
\]

It reuses the validated 3D scalar/tensor diffusion assembly and supports scalar or
symmetric positive-definite tensor permittivity, spatial and additive per-cell
free-charge density, fixed nodal or exterior-face potentials, and prescribed
outward displacement D dot n. Pure-Neumann displacement data use the existing
component compatibility check and zero-volume-mean gauge. Unassigned exterior
faces have zero normal displacement. The result includes nodal potential, per-cell E and D, and
stored energy U = 1/2 integral(E dot D) dV.

~~~sh
./build/nova_electrostatics_3d
~~~

An experimental [fixed-electrode charge diagnostic](docs/ELECTROSTATICS_3D.md#fixed-electrode-charge-diagnostic-experimental)
integrates element-averaged displacement on validated prescribed-potential
boundary faces. It reports dielectric-outward flux and oppositely signed
conductor charge, but is not an equilibrated Gauss-law flux.

An experimental [multi-terminal capacitance matrix extraction](docs/ELECTROSTATICS_3D.md#multi-terminal-capacitance-matrix-experimental)
performs unit-voltage basis solves and reports the row-major Maxwell matrix
plus its reciprocity defect. It uses the same approximate element-averaged
boundary charge diagnostic, not discrete electrode reactions.

The analytical suite covers a parallel-plate capacitor (C = 2 F, U = 100 J),
an anisotropic affine field (U = 36.5 J), and an exact one-free-DOF tetrahedral
charge benchmark. The installed SDK consumer also verifies the capacitor solution, fixed-electrode charge API and multi-terminal capacitance matrix.
This is an experimental electrostatic Poisson interface, not a general Maxwell
solver; floating conductors, surface-charge unknowns, AC/transient Maxwell,
and equilibrated displacement recovery are not implemented. See
[3D electrostatics formulation and limits](docs/ELECTROSTATICS_3D.md).

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

Examples: examples/transient_heat_2d/main.cpp (P2 cooling plate), examples/transient_flux_2d/main.cpp (P2 Crank–Nicolson boundary exchange), and examples/adaptive_heat_2d/main.cpp (error-controlled P2 Crank–Nicolson heating). See [adaptive 2D API, caching behavior and limitations](docs/ADAPTIVE_HEAT_2D.md).

Limitations: no temperature-dependent material coefficients, adaptive mesh refinement, time integrators above order two, or independently reconstructed Dirichlet normal heat flux. Adaptive stepping reuses the spatial operators and dt-dependent base systems, but time-dependent Robin matrices and constrained systems are still assembled at their required time levels. This prototype is not validated for engineering design use.

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

### Verified 3D L2 and H1 spatial-error measurement

The P1 tetrahedral verification utility integrates exact-solution L2 and H1
errors with 125-point Duffy-mapped Gaussian quadrature. Automated analytical
patch, orientation, quadratic interpolation and manufactured Poisson tests
check the expected spatial convergence rates.

See [3D error-norm formulation and verification](docs/VERIFICATION_3D_ERROR_NORMS.md).

### Linear and nonlinear solver capabilities
- CSR sparse matrix-vector products.
- Conjugate Gradient for symmetric positive-definite systems.
- BiCGSTAB and restarted GMRES for general nonsymmetric systems.
- Optional Jacobi diagonal preconditioning for all three iterative linear solvers.
- Natural-order zero-fill incomplete Cholesky (IC(0)) preconditioning for Conjugate Gradient on symmetric systems with resolved positive pivots.
- Newton nonlinear solver with GMRES or BiCGSTAB linear subsolves.
- Armijo backtracking, including rejection of non-finite trial residuals.
- Explicit breakdown and convergence reporting.
- True-residual recomputation before linear convergence is accepted.

Select IC(0) for a CG solve with
`options.preconditioner = nova::linalg::PreconditionerType::IncompleteCholesky0`.
IC(0) is an unshifted natural-order zero-fill factor: it rejects asymmetric
matrices and unresolved/non-positive pivots rather than silently changing the
system or falling back. Some SPD matrices do not admit a usable natural-order
IC(0) factor; use Jacobi or no preconditioner for those cases. IC(0) is CG-only;
GMRES and BiCGSTAB continue to accept None or Jacobi. Use an instance of `nova::linalg::CachedConjugateGradient` when solving
multiple right-hand sides. It retains one IC(0) factor and reuses it only when
matrix dimensions, sparsity and every stored coefficient are exactly
unchanged; any coefficient or pattern change triggers a rebuild. The fixed-step
2D and 3D transient solvers own one cache per integration. Time-dependent Robin
coefficients therefore refresh the factor when they change the system matrix.
Cache instances are stateful and not thread-safe. Adaptive 2D and 3D step
doubling use separate full-step and half-step cache instances across trials.
The 2D solver also reevaluates time-dependent loads and Robin matrices at each
time level; exact constrained-CSR identity decides whether IC(0) can be reused.

The 2D P1/P2 stationary solver supports compatible pure-Neumann components with a zero-area-mean gauge; it remains scalar and isotropic. The 3D stationary solver also supports balanced floating components and spatial SPD conductivity tensors. Multilevel preconditioning, error-driven 3D mesh adaptation and a general multiphysics field registry remain outside this mainline configuration. The nonlinear 1D example is a verification model, not yet a general nonlinear material framework.

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
./build/nova_nonlinear_diffusion_1d
./build/nova_electrostatics_3d
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
