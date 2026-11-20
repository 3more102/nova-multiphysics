# Transient 3D heat conduction (experimental C++ API)

## Governing equation

NOVA now has a time-dependent 3D **P1 conforming tetrahedral** heat solver:

```text
C(x) dT/dt - div(K(x) grad T) = Q(x)
```

- `C = rho * c_p > 0`: volumetric heat capacity, J/(m³ K).
- `K`: scalar positive conductivity or symmetric-positive-definite tensor, W/(m K).
- `Q`: signed volumetric generation, W/m³.
- `T`: temperature (kelvin when inputs use SI).
- Outward Neumann heat flux is positive for energy leaving the solid.
- Robin: outward flux = `h * (T - T_ambient)`.

The spatial operators are **consistent** P1 mass and conductivity stiffness matrices. Both are assembled with the existing positive degree-two tetrahedral quadrature; exterior Neumann/Robin face terms use the existing degree-two triangular quadrature. The code reuses `DiffusionMaterials3D`, `integrate_p1_volume_3d`, mesh validation and exclusive exterior-face assignment from stationary conduction.

## Time discretization

Given the assembled heat-capacity matrix `M`, thermal stiffness and convection matrix `K`, and total fixed load `F`:

```text
[M + theta * dt * K] T_(n+1) =
  [M - (1-theta) * dt * K] T_n + dt * F

theta = 1   Backward Euler
theta = 1/2 Crank–Nicolson
```

`F` includes volume heating, **negative** outward Neumann flux, and the Robin ambient load. Initial nodal values are projected to the fixed Dirichlet constraints. Each step uses symmetric elimination with identity rows on constrained nodes. The sparse constrained matrix is constructed **once** and the preceding solution warm-starts CG. Positive `C` makes the time-step system nonsingular even for a completely insulated mesh.

The per-step discrete energy budget uses `M * ones` as the heat-content integration weights. Dirichlet input is calculated from the unconstrained system residual on prescribed nodes, **not** from raw P1 face flux. For constant coefficients and loads,

```text
Delta(Integral C*T dV) =
    source_energy + neumann_energy + robin_energy + dirichlet_energy
    + iterative_balance_residual
```

## Minimal API usage

```cpp
#include "nova/fem/transient_heat_3d.hpp"

using namespace nova::fem;
TransientHeat3DProblem problem;
problem.diffusion.mesh = make_box_mesh_3d(0,1,0,1,0,1,2,2,2).mesh;
problem.diffusion.cell_diffusivity.assign(
    problem.diffusion.mesh.elements.size(), 10.0);
problem.diffusion.source = [](double, double, double) { return 1000.0; };
problem.volumetric_capacity = [](double, double, double) { return 3.0e6; };
problem.initial_temperature.assign(problem.diffusion.mesh.nodes.size(), 300.0);

TransientHeat3DOptions options;
options.time_step = 0.1;
options.steps = 10;
options.scheme = TransientHeat3DScheme::BackwardEuler;
options.solver_options.preconditioner =
    nova::linalg::PreconditionerType::Jacobi;

const auto output = solve_transient_heat_3d(problem, options);
```

The `diffusion` member also accepts the existing steady 3D `dirichlet_constraints`, `dirichlet_faces`, `neumann_faces`, `robin_faces`, `cell_source`, scalar and tensor conductivity representations. Heat capacity must be provided **either** as `volumetric_capacity` **or** a `cell_volumetric_capacity` vector; both simultaneously are invalid.


### Zero and cell-only source terms

The optional `diffusion.source` function defaults to zero spatial generation
when absent. A `diffusion.cell_source` vector alone can provide region-specific
heating, without a dummy zero-valued lambda; both fields absent describe
unforced conduction. Finite-value checks still apply at every quadrature point.

### Analytical temporal-order verification

A unit tetrahedron with three fixed vertices leaves one free temperature DOF. For the reference P1 basis, the free-node entries are M33 = 1/60 and K33 = 1/6, so the semidiscrete equation is dT/dt = -10 T. Starting from T(3) = 10, the exact value at t = 3.2 is 10 exp(-2). The regression compares time steps of 0.05 and 0.025: Backward Euler's global error ratio approaches 2, while Crank–Nicolson's approaches 4. It checks every accepted step's discrete energy residual and confirms the second-order result is more accurate.

## Reproduce and tests

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

`tests/test_transient_heat_3d.cpp` exercises:

- Uniform insulated heating: `C=4`, `Q=8`, `T(0)=300` gives `T(t)=300+2t`, exactly represented by P1, for both schemes.
- Three fixed vertices of a unit tetrahedron: the remaining DOF has `M_33=1/60`, `K_33=1/6`. Starting at 10 with `dt=0.1`, the one-step exact algebraic results are 5 (Backward Euler) and 10/3 (Crank–Nicolson).
- Neumann-flux heat loss, Robin convection, mixed tensor conductivity with a spatial heat capacity, and numerical energy-budget closure.
- Heterogeneous cellwise SPD tensor conductivity with cellwise capacity and cell-only source, plus zero-source insulated equilibrium.
- Multi-tetrahedron heterogeneous cellwise SPD tensors with uniform heating, exact heat-content growth, and invalid tensor/cell-array rejection.
- Analytical temporal refinement verifies first-order Backward Euler and second-order Crank–Nicolson at a nonzero physical start time.
- Rejection of nonpositive capacities, bad source/initial fields, duplicate boundary conditions and invalid time-step options.

These checks verify the mathematical discretization but **do not** constitute experimental or industrial validation.

## Streaming snapshots and ParaView export

The `TransientHeat3DOptions::on_step` observer receives
`(step_index, physical_time, nodal_temperatures)` after initial Dirichlet
projection at step **0**, then after each converged step. Its field reference
is valid only during the callback. Observer exceptions propagate without
silently restarting or skipping a step.

By default, `result.linear_solvers[i].solution` still contains each full
nodal state for backward compatibility. Setting
`options.retain_step_solutions = false` keeps CG convergence diagnostics
but clears those stored vectors; the final `result.temperature` is still
returned. This avoids an O(number of steps × number of nodes) full-field
result history, although on-disk VTK frames themselves consume storage.

```sh
./build/nova_transient_heat_3d --vtk-series thermal_case.pvd
```

This writes `thermal_case.pvd` and `thermal_case_000000.vtk` through
`thermal_case_000005.vtk`, with times 0, 20, 40, 60, 80 and 100 seconds.
No extra solve is run per frame. The PVD writer requires unique, safe ASCII
`.vtk` basenames and strictly increasing finite physical times; invalid
collections are rejected before XML output. The collection file is
written only after every VTK frame and step has completed successfully,
but previously written frames are not removed if later export fails.

`nova_transient_heat_3d_pvd_verify` parses all generated XML and VTK frames
with Python's standard library and compares the final field range to a
separate monolithic solve. Unit tests cover the observer's initial and final
fields, Dirichlet projection, error propagation and optional history elision.

## Explicit limits and next milestone

Initial 3D support uses **time-independent** conductivity, capacity, forcing and boundary functions; a constant time step; conforming straight-sided P1 tetrahedral meshes; and an in-memory C++ API. An experimental, separate fixed-mesh 3D adaptive timestep driver is documented in [Adaptive transient 3D heat](ADAPTIVE_HEAT_3D.md), but temperature-dependent properties, latent heat, radiation and moving meshes remain unsupported. The default result history retains every CG solution vector for API compatibility; the optional streaming mode elides them. Next: reuse assembled operators for adaptive 3D trials, spatial error-driven refinement, and independent temporal-convergence benchmarks on spatially varying manufactured solutions.
