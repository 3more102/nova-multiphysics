# Nonlinear transient thermal prototype: frozen conductivity

## Governing equation

The first nonlinear 2D thermal prototype solves, in weak Galerkin form,

`rho*c(x,y) * dT/dt - div(k(x,y,T) grad T) = Q(x,y,t)`

on fixed P1/P2 triangular H1 meshes. Positive capacity depends on position
but not temperature. Conductivity is an isotropic, strictly positive
user-defined callback `k(x,y,T)`.

This is a **Backward Euler + Picard** reference implementation, not a Newton
method. Crank–Nicolson and adaptive time stepping remain available only in
the linear thermal path. A temperature-dependent mass/capacity model is not
yet implemented.

## Solution procedure

For each time step:

1. Start from the previous accepted H1 DOF state.
2. Interpolate the current Picard iterate inside each triangle, then freeze
   conductivity at `k(x,y,T_iter(x,y))` at all thermal quadrature points.
3. Call the existing symmetric P1/P2 Backward Euler linear thermal solver
   with that conductivity and the previous *accepted* state as initial data.
   Volume sources, P1/P2 Dirichlet temperatures, Neumann fluxes and Robin
   convection use the existing boundary-condition semantics.
4. Compare the linear solve result with the Picard iterate using a
   maximum scaled nodal update:

   `max_i |T_trial[i]-T_iter[i]| /
     (atol + rtol * max(|T_trial[i]|, |T_iter[i]|)) <= 1`.

5. Accept the trial solution if converged; otherwise update
   `T_iter += omega * (T_trial-T_iter)` with `0<omega<=1`.
6. Throw if the configured maximum iteration budget is exhausted. No
   silently accepted partially converged time step is returned.

`PicardHeat2DResult` stores the final H1 DOFs, physical end time,
iterations/scaled updates per accepted time step, and all inner linear
solver reports.

## C++ example

```cpp
nova::fem::NonlinearHeat2DProblem problem;
problem.thermal.mesh = mesh;
problem.thermal.volumetric_capacity =
    [](double, double) { return 1.0; };
problem.thermal.source =
    [](double, double, double) { return 0.0; };
problem.thermal.initial_temperature =
    [](double x, double y) { return 1.0 + x * x + y * y; };
problem.conductivity =
    [](double, double, double temperature) {
        return 1.0 + temperature * temperature;
    };

nova::fem::PicardHeat2DOptions options;
options.order = nova::fem::H1TriangleOrder::P2;
options.time_step = 0.05;
options.steps = 4;
const auto result = nova::fem::solve_picard_heat_2d(problem, options);
```

The `thermal.conductivity` callback, if supplied, is not used: it is
replaced at each Picard iteration by the temperature-dependent callback.

## Verification coverage

Automated tests include P1/P2 constant-state preservation, an affine
space-time manufactured solution with `k=1+0.2*T`, an evolving nonuniform
P2 field, nonpositive conductivity rejection, control validation and
iteration-budget failure. These are regression and manufactured-solution
checks, not independent industrial validation.

## Known limits and priorities

- The conductivity lookup currently scans all mesh triangles for each
  assembly query, so cost increases strongly with element count. Replace
  it with element-aware quadrature coefficient evaluation before scale-up.
- No consistent nonlinear tangent/Jacobian, globalization strategy beyond
  simple Picard under-relaxation, or guaranteed convergence for strong
  nonlinearity.
- No nonlinear Crank–Nicolson/adaptive controller and no independent
  nonlinear energy-balance or exact-flux certification.
- Scalar isotropic conductivity only; no temperature-dependent capacity,
  phase change, latent heat, anisotropic materials, or coupled multiphysics.
- Inherited iterative solver tolerances and scaled DOF changes are
  numerical heuristics, not engineering error certificates.

Next: move coefficient evaluation into element assembly, add residual-based
convergence tests, and implement consistent-tangent Newton for `k(T)`.
