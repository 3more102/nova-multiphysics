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

## Streaming snapshots and memory control

An optional `TransientHeat3DOptions::on_state` callback receives the
**Dirichlet-projected initial state** at step 0, then the completed temperature
field at each step, together with the physical time. Its vector reference is
valid only during the callback. This supports external VTK/PVD frame writers
without accumulating all snapshots in memory.

By default, `result.linear_solvers[n].solution` contains the physical
temperature at step `n+1`, preserving the original API. Set
`options.retain_step_solutions = false` for long runs: all per-step solver
metadata and diagnostics remain available, but each stored solver's
`solution` is empty; `result.temperature` still contains the final field.
This reduces retained temperature history from O(steps × nodes) to O(nodes),
excluding user-managed copies in the callback. The diagnostics and solver
metadata themselves still consume O(steps) space.

```cpp
options.retain_step_solutions = false;
options.on_state = [](std::size_t step, double time,
                      const std::vector<double>& temperature) {
    // Stream or serialize this temperature field before returning.
};
```

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
- Rejection of nonpositive capacities, bad source/initial fields, duplicate boundary conditions and invalid time-step options.
- Streaming callback receives projected state at t=0 and every completed step without retaining solver solution vectors.
- Semidiscrete one-DOF exponential decay confirms first-order Backward Euler and second-order Crank–Nicolson time convergence under time-step refinement.

These checks verify the mathematical discretization but **do not** constitute experimental or industrial validation.

## Explicit limits and next milestone

Initial 3D support uses **time-independent** conductivity, capacity, forcing and boundary functions; a constant time step; conforming straight-sided P1 tetrahedral meshes; and an in-memory C++ API. No 3D adaptive stepping, temperature-dependent properties, latent heat, radiation, moving meshes, transient Gmsh CLI or ParaView time-series driver is included in this PR. Streaming snapshots and disabling retained solver solutions are now supported. Next: a Gmsh-tagged transient CLI with `.pvd` export, plus temporal-convergence verification on a spatially varying manufactured solution with multiple free degrees of freedom.
