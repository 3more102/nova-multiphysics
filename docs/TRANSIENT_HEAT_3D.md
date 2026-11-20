# Transient 3D heat conduction (experimental C++ API)

## Governing equation

NOVA now has a time-dependent 3D **P1 conforming tetrahedral** heat solver:

```text
C(x) dT/dt - div(K(x) grad T) = Q(x,t)   [static capacity]
d(C(x,t) T)/dt - div(K(x,y,z,t) grad T) = Q(x,t)   [dynamic capacity]
```

- `C = rho * c_p > 0`: volumetric heat capacity, J/(m³ K).
- `K`: scalar positive conductivity or symmetric-positive-definite tensor, W/(m K).
- `Q`: signed volumetric generation, W/m³.
- `T`: temperature (kelvin when inputs use SI).
- Outward Neumann heat flux is positive for energy leaving the solid.
- Robin: outward flux = `h * (T - T_ambient)`.

The spatial operators are **consistent** P1 mass and conductivity stiffness matrices. Volume terms use the existing positive degree-two tetrahedral quadrature; exterior Neumann/Robin face terms use degree-two triangular quadrature. Static scalar/tensor material validation shares `DiffusionMaterials3D`; dynamic conductivity callbacks are evaluated at every material quadrature point and time level. Mesh validation and exclusive exterior-face assignment are shared with stationary conduction.

## Time discretization

Given the assembled heat-capacity matrix `M`, the time-level thermal operator `K(t)` (conductivity plus static and time-dependent Robin convection), and total load `F(t) = F_static + F_time(t)`:

```text
[M(t_(n+1)) + theta * dt * K(t_(n+1))] T_(n+1) =
  M(t_n) T_n - (1-theta) * dt * K(t_n) T_n
  + dt * [(1-theta) * F(t_n) + theta * F(t_(n+1))]

theta = 1   Backward Euler
theta = 1/2 Crank–Nicolson
```

`F` includes volume heating, **negative** outward Neumann flux, and Robin ambient loads. The optional `problem.time_dependent_source(x,y,z,t)` adds a time-varying volumetric source to static spatial and cell sources. Backward Euler samples loads at the new endpoint; Crank–Nicolson averages endpoint loads. Time-dependent scalar/tensor conductivity and volumetric-capacity callbacks are also sampled at both material time levels. When capacity varies, the solver uses conservative enthalpy storage `d(C*T)/dt`, so the discrete heat content is `1^T M(t) T(t)`; it does not silently reinterpret the callback as `C(t) dT/dt`. Initial nodal values are projected to the fixed Dirichlet constraints. Each step uses symmetric elimination with identity rows on constrained nodes. Exact endpoint-matrix comparison controls operator and IC(0) reuse. Positive capacity makes the time-step system nonsingular even for a completely insulated mesh.

The per-step discrete energy budget uses `M * ones` as the heat-content integration weights. Dirichlet input is calculated from the unconstrained system residual on prescribed nodes, **not** from raw P1 face flux. For constant coefficients and loads,

```text
Delta(Integral C(t)*T dV) =
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

The `diffusion` member also accepts the existing steady 3D `dirichlet_constraints`, `dirichlet_faces`, `neumann_faces`, `robin_faces`, `cell_source`, scalar and tensor conductivity representations. Static heat capacity must be provided **either** as `volumetric_capacity` **or** a `cell_volumetric_capacity` vector. Alternatively, use `time_dependent_volumetric_capacity(x,y,z,t)`; it is mutually exclusive with both static forms.


### Zero and cell-only source terms

The optional `diffusion.source` function defaults to zero spatial generation
when absent. A `diffusion.cell_source` vector alone can provide region-specific
heating, without a dummy zero-valued lambda; both fields absent describe
unforced conduction. Finite-value checks still apply at every quadrature point.

### Analytical temporal-order verification

A unit tetrahedron with three fixed vertices leaves one free temperature DOF. For the reference P1 basis, the free-node entries are M33 = 1/60 and K33 = 1/6, so the semidiscrete equation is dT/dt = -10 T. Starting from T(3) = 10, the exact value at t = 3.2 is 10 exp(-2). The regression compares time steps of 0.05 and 0.025: Backward Euler's global error ratio approaches 2, while Crank–Nicolson's approaches 4. It checks every accepted step's discrete energy residual and confirms the second-order result is more accurate.

### Time-dependent Dirichlet, Neumann and Robin data

The transient API accepts time callbacks `f(x,y,z,t)` for nodal
Dirichlet constraints, exterior Dirichlet faces, exterior Neumann face fluxes,
and exterior Robin transfer/ambient fields. Dirichlet values are imposed at
the initial physical time and at each step endpoint. Neumann loads are
integrated over triangular faces using the existing degree-two face
quadrature.

Backward Euler samples time-dependent Neumann flux and Robin data at the new
endpoint. Crank–Nicolson evaluates both endpoints and uses trapezoidal
integration. For Robin convection, both the transfer matrix
`R(t) = integral_face h(x,t) phi_i phi_j dA` and ambient load
`g(t) = integral_face h(x,t) T_ambient(x,t) phi_i dA` are evaluated at the
required time levels. The API field is `problem.time_dependent_robin_faces`:

```cpp
problem.time_dependent_robin_faces.push_back({
    face,
    [](double, double, double, double time) { return 12.0 + 0.1 * time; },
    [](double, double, double, double time) { return 290.0 + 2.0 * time; }
});
```

Time-dependent Dirichlet values change the prescribed vector but not the
constrained matrix pattern. Time-dependent Neumann data and Robin ambient
temperature change only the right-hand side. A changed Robin transfer
coefficient changes the time-level matrix. Time-dependent scalar/tensor
conductivity and volumetric capacity callbacks are evaluated at both endpoints:
Backward Euler uses the old capacity matrix on the old state and the new
capacity/conductivity operators implicitly; Crank–Nicolson uses both endpoint
material matrices. A capacity callback represents conservative
`d(C*T)/dt) storage. Exact constrained-CSR comparison controls operator and
IC(0) reuse. The same endpoint evaluation is used by adaptive full-step and
half-step trials.

### Time-dependent material callbacks

Choose one dynamic conductivity representation:
`time_dependent_diffusivity(x,y,z,t)` for a positive scalar field, or
`time_dependent_diffusivity_tensor(x,y,z,t)` for an SPD tensor. Each is
mutually exclusive with all static scalar, tensor and cellwise conductivity
fields. Dynamic capacity is supplied separately:

```cpp
problem.time_dependent_diffusivity =
    [](double, double, double, double time) { return 10.0 * (1.0 + 0.01 * time); };
problem.time_dependent_volumetric_capacity =
    [](double, double, double, double time) { return 3.0e6 * (1.0 + 0.005 * time); };
```

The capacity callback is mutually exclusive with both static capacity forms.
For a changing capacity, NOVA discretizes the conservative storage term
`d(C(x,t)T)/dt); this makes the reported heat-content change consistent with
the discrete source and boundary energy budget. Temperature-dependent
constitutive laws (callbacks depending on the unknown temperature) are not
implemented.

## Reproduce and tests

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

`tests/test_transient_heat_3d.cpp` exercises:

- Uniform insulated heating: `C=4`, `Q=8`, `T(0)=300` gives `T(t)=300+2t`, exactly represented by P1, for both schemes.
- A uniform ramp source `Q_time(t)=t` verifies Backward Euler endpoint sampling, Crank–Nicolson trapezoidal integration and source-energy closure.
- Three fixed vertices of a unit tetrahedron: the remaining DOF has `M_33=1/60`, `K_33=1/6`. Starting at 10 with `dt=0.1`, the one-step exact algebraic results are 5 (Backward Euler) and 10/3 (Crank–Nicolson).
- Neumann-flux heat loss, static and time-dependent Robin convection, mixed tensor conductivity with a spatial heat capacity, and numerical energy-budget closure.
- Closed-form one-free-DOF Backward Euler and Crank–Nicolson recurrences for time-varying Robin transfer/ambient data, including per-step Robin energy.
- Skewed, non-axis-aligned tetrahedron checks of time-dependent Neumann and Robin face quadrature against independent one-free-DOF recurrences and per-step energy closure.
- Adaptive full-step/two-half-step integration with time-dependent Robin data and energy-budget closure.
- Heterogeneous cellwise SPD tensor conductivity with cellwise capacity and cell-only source, plus zero-source insulated equilibrium.
- Multi-tetrahedron heterogeneous cellwise SPD tensors with uniform heating, exact heat-content growth, and invalid tensor/cell-array rejection.
- Analytical temporal refinement verifies first-order Backward Euler and second-order Crank–Nicolson at a nonzero physical start time.
- Time-dependent scalar and SPD-tensor conductivity plus time-dependent capacity, with conservative heat-content diagnostics and BE/CN convergence-order checks.
- Rejection of nonpositive static/dynamic capacities, invalid scalar/tensor dynamic conductivity, mixed static/dynamic material representations, bad source/initial fields and invalid time-step options.

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

Initial 3D support includes static or prescribed time-varying scalar/tensor conductivity and volumetric capacity, additive time-dependent volumetric sources, Dirichlet constraints/faces, Neumann face fluxes and Robin transfer/ambient functions; fixed-step and experimental fixed-mesh adaptive time integration; conforming straight-sided P1 tetrahedral meshes; and an in-memory C++ API. Dynamic capacity uses conservative enthalpy storage. Temperature-dependent constitutive laws, latent heat, radiation, moving meshes and spatial error-driven tetrahedral refinement remain unsupported. The default result history retains every CG solution vector for API compatibility; optional streaming elides them. Next: spatial error-driven refinement and independent temporal-convergence benchmarks on spatially varying manufactured solutions.
