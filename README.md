# NOVA Multiphysics

**Engineering the Physics of Tomorrow.**

NOVA is an independent C++20 multiphysics simulation platform under development. Its long-term objective is a commercial, general-purpose CAE environment for geometry, meshing, finite-element modeling, numerical solution, multiphysics coupling, and scientific visualization.

> **Development status:** Early engineering prototype. It is not production-certified and must not be used as the sole basis for engineering design decisions.

[Project status and PR dependencies](docs/PROJECT_STATUS.md) · [Engineering decisions](docs/DECISIONS.md) · [Release criteria](docs/RELEASE_CRITERIA.md)

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
- Pure-Neumann loads checked per connected component, with a zero area-mean P1 or P2 field. Positive Robin transfer also anchors a component.
- L2 and H1 error norms with Duffy-transformed Gaussian quadrature.
- Orientation-independent triangle handling and scale-aware degenerate-element rejection.

### 3D scalar diffusion
The tetrahedral P1 H1 kernel solves

\[
-\nabla\cdot(k\nabla u)=f \quad \text{in }\Omega
\]

on validated four-node tetrahedral meshes using:
- Positive four-point degree-two tetrahedron quadrature.
- Three-point degree-two triangular-face quadrature.
- Signed-Jacobian physical P1 gradients and orientation-independent volume integration.
- Spatially varying scalar or symmetric positive-definite tensor diffusivity and volumetric sources.
- Strong Dirichlet conditions on vertices or complete exterior faces.
- Neumann outward-flux integration and Robin exchange on validated exterior faces.
- Homogeneous natural (zero-flux) conditions on unassigned faces.
- Symmetric Dirichlet elimination and conjugate-gradient solution.
- Pure-Neumann components with compatible source/flux loads and zero volume-weighted mean.
- Affine harmonic patch, variable-coefficient, mixed-boundary and nodal-convergence tests.

Example: `examples/diffusion_3d/main.cpp`.

Pure-Neumann compatibility is checked separately on each connected component without Dirichlet or positive Robin anchoring. The returned field on each such component is normalized to zero volume-weighted mean.

For direction-dependent materials, set `diffusivity_tensor` on a stationary problem or `conductivity_tensor` on a transient problem. Supply the six components in global Cartesian coordinates and leave the scalar callback empty. Fixed and adaptive time stepping use the same material assembly. [Tensor conventions and supported range](docs/NUMERICAL_CONVENTIONS.md#tetrahedral-conductivity-tensors)

### Import external tetrahedral meshes (Gmsh MSH 2.2 ASCII)

Use `nova::io::read_gmsh22_tetrahedra(std::istream&)` from
`nova/io/gmsh22_tetrahedral.hpp` to import a mesh generated externally
by Gmsh. The reader accepts linear four-node tetrahedra (Gmsh type 4),
optional surface triangles (type 2), and line/point records. It remaps
arbitrary positive Gmsh node IDs to zero-based NOVA P1 vertex indices,
preserves per-tetrahedron and per-exterior-face physical IDs, and
reorients boundary faces outward. Unassigned exterior faces receive
physical ID 0. The returned `mesh` can be passed directly to
`ScalarDiffusion3DProblem` or the transient 3D thermal problem.

Only **MSH 2.2 ASCII, double precision, first-order tetrahedra** are
supported. MSH 4.x, binary meshes, higher-order cells and other volume
cell types are rejected instead of silently discarded. Topology,
coordinates, element references and exported surface labels are checked;
a mesh with isolated vertices, nonmanifold faces, degeneracy or
surface triangles on interior faces is rejected.

For example, generate a first-order 3D mesh in Gmsh using
`gmsh geometry.geo -3 -format msh2 -o model.msh`. Then:

~~~cpp
std::ifstream mesh_file("model.msh");
auto imported = nova::io::read_gmsh22_tetrahedra(mesh_file);
// imported.mesh; imported.exterior_faces;
// imported.exterior_physical_tags; imported.volume_physical_tags
~~~

The caller must explicitly assign boundary conditions and material
properties from the physical IDs; the importer does not infer them.

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
physical boundary groups through a general registry. The H1 transfer and remeshing helpers preserve P1/P2 fields on supported refinements.

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

Limitations: no temperature-dependent coefficients, multi-cycle adaptive mesh refinement, time integrators above order two, or independently reconstructed Dirichlet normal heat flux. The 2D adaptive prototype reassembles FEM matrices for each trial; the new 3D adaptive path reuses static mass/stiffness operators but still rebuilds time-step-dependent systems and time-dependent loads. This prototype is not validated for engineering design use.

### Transient 3D heat

The tetrahedral P1 solver supports consistent thermal mass, fixed-step Backward Euler and Crank–Nicolson, spatially varying positive capacity and conductivity, time-dependent Dirichlet, outward Neumann and Robin face conditions, exact nodal restarts and per-step energy accounting. The stationary and transient solvers share quadrature and exterior-face validation. Insulated transient components preserve their initial heat content.

The adaptive 3D path adds full-step/two-half-step control for Backward Euler and Crank–Nicolson, automatic acceptance/rejection, RMS absolute/relative tolerances, exact P1 restart states, combined energy diagnostics, and one-time assembly of static mass/stiffness operators per study. For a tensor-conductivity cube, independent scikit-fem assembly and SciPy sparse solves check each trial. SciPy Radau checks the terminal temperatures on a different time grid. Run `nova_anisotropic_heat_3d cn 4 1e-6` to reproduce the tighter Crank–Nicolson case. [Adaptive verification and timings](docs/THERMAL_3D_VERIFICATION.md#adaptive-tensor-heat)

Run `nova_transient_heat_3d 8 20` for a cooling cube with a known analytical solution, or `nova_adaptive_heat_3d` for adaptive P1 heating. The optional third argument `--csv` writes all final nodal temperatures. [Fixed-step verification results and reproduction commands](docs/THERMAL_3D_VERIFICATION.md)

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
- Zero-fill incomplete Cholesky for CG on symmetric systems with suitable pivots.
- Newton nonlinear solver with GMRES or BiCGSTAB linear subsolves.
- Armijo backtracking, including rejection of non-finite trial residuals.
- Explicit breakdown and convergence reporting.
- True-residual recomputation before linear convergence is accepted.

Use `solver_options.preconditioner = nova::linalg::PreconditionerType::IncompleteCholesky0` to select IC(0). Use `nova_diffusion_3d ic0 24` or `nova_diffusion_3d jacobi 24` to compare the same manufactured problem. [Iteration counts, true residuals and full-process timings](docs/PROJECT_STATUS.md#19-november-2026--zero-fill-incomplete-cholesky)

Open milestones include factor reuse across heat steps, P2 tetrahedra, adaptive 3D meshing and a general multiphysics field registry. A general nonlinear material framework would extend the existing 1D verification model.

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
./build/nova_diffusion_3d
./build/nova_transient_heat_2d
./build/nova_transient_heat_3d 8 20
./build/nova_transient_flux_2d
./build/nova_adaptive_heat_2d
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
