# Experimental implicit 3D transient electrothermal coupling

Public API: `nova/physics/transient_electrothermal_3d.hpp`,
`solve_transient_electrothermal_3d`.

## Physics and units

For each time step `t_n -> t_(n+1)`, this segregated algorithm
iterates the *same* tetrahedral mesh for quasi-static DC electrical conduction,
P1 Joule dissipation, and Backward Euler thermal diffusion:

```text
-div(sigma(T_trial) grad V) = impressed_current_density  [A/m^3]
J = -sigma(T_trial) grad V                             [A/m^2]
q_J = E dot J                                          [W/m^3]
C (T_(n+1)-T_n)/dt - div(k(T_trial) grad T_(n+1))
    = q_J + q_base                                     [W/m^3]
```

`C` is the original thermal solver's positive volumetric heat capacity
in J/(m^3 K). The standard transient FEM solver also supports conservative
time-dependent `C(t)` and time-dependent boundary/load data; these are not
dependent on temperature. This API deliberately **only** uses Backward Euler,
not Crank-Nicolson. Electrical potential, E, current density and q_J are
recomputed on every nonlinear Picard iteration. No displacement current.

A caller-supplied optional electrical/thermal material law scales the
reference **per-cell scalar or SPD-tensor** coefficient by a positive
dimensionless factor evaluated at the four-node mean trial temperature.
Without the law, the original static electrical/thermal material is used.
Each scaling law is finite and strictly positive. A thermal k(T) law cannot
be combined with the heat solver's time-dependent conductivity callback.

`diffusion.source`, `diffusion.cell_source`, `time_dependent_source`
and the computed Joule density are additive; no caller-owned data is altered.
Neumann and Robin heat conditions retain their sign conventions. The returned
`heat_budget` is the discrete transient weak-form energy budget in joules,
not an independently certified continuum boundary flux.

## Numerical algorithm

For each step and each Picard iteration:

1. Initialize `T_trial` from the last accepted temperature `T_n`.
2. Update electrical `sigma(T_trial)`; solve the stationary DC field and
   compute cellwise `q_J`.
3. Update thermal `k(T_trial)`, add `q_J` to the original thermal cell
   sources, and solve exactly one *implicit* Backward Euler heat step
   **from the last accepted** `T_n`. Do not advance the thermal initial
   state during failed nonlinear trial iterations.
4. Check the **unrelaxed** fixed-point defect
   `|T_solved - T_trial| / (abs_tol_K + rel_tol*max(|T_solved|,|T_trial|))`.
5. If it is >1, relax `T_trial` toward `T_solved` and retry. After an
   unsuccessful iteration cap, return `converged=false` with
   `failed_step` and the last accepted thermal state; never commit a
   partly converged step.

`accepted_steps` contains end times, nonlinear iteration count,
raw and normalized residuals, Joule watts and the transient heat budget.
`last_electrical` is the electrical solve from the **last accepted**
temperature fixed-point iteration. It is within the configured numerical
fixed-point tolerance, not evaluated at an independently re-solved exact
final temperature.

## Analytic checks

For a unit cube with insulating thermal boundaries, V(0)=0, V(1)=10 V,
reference sigma=2 S/m and C=2 J/(m^3 K), the Joule source is 200 W/m^3.
With T(0)=300 K and constant sigma, `T(t)=300+100t` exactly under
the P1 backward-Euler spatial discretization. Each 0.02 s step generates
4 J of Joule energy in the unit cube.

With `sigma(T)=2*[1+0.01*(T-300)]` S/m, uniform insulated temperature,
and the same capacity and voltage, Backward Euler obeys the independent
closed-form scalar relation, with `D=T-300` and dt=0.05 s:

```text
D_(n+1) = (D_n + 100*dt) / (1 - dt).
```

Because temperature stays uniform, a positive cellwise k(T) changes the
thermal stiffness but does not alter this exact uniform-temperature
reference. The tests also check tensor conductivity, additive heat sources,
input rejection and iteration-cap behavior.

## Example C++ usage

```cpp
#include "nova/physics/transient_electrothermal_3d.hpp"

nova::physics::TransientElectrothermal3DProblem problem;
problem.electrical = electrical_dc_problem;
problem.thermal = transient_heat_problem;  // same mesh and ordering
problem.electrical_conductivity_scale = [](double t) {
    return 1.0 + 0.003 * (t - 300.0);
};
problem.thermal_conductivity_scale = [](double t) {
    return 1.0 + 0.01 * (t - 300.0);
};
nova::physics::TransientElectrothermal3DOptions options;
options.thermal_time.time_step = 0.01;
options.thermal_time.steps = 100;
options.max_coupling_iterations = 40;
const auto result = nova::physics::solve_transient_electrothermal_3d(problem, options);
if (!result.converged) {
    // result.failed_step and result.failed_max_normalized_residual
}
```

## Explicit limitations and next steps

- Fixed time step and Backward Euler only; no step-doubling, time-adaptive
  accept/reject, nonlinear Newton Jacobian or adaptive relaxation.
- The solver solves a new DC and thermal FEM system at every Picard trial.
  It does not cache or reuse all operators across steps; large simulations
  may be expensive.
- Temperature-dependent scalar/tensor conductivity uses **cell-average P1
  nodal temperature**, not temperature evaluated at quadrature points.
- No Seebeck/Peltier coupling, thermal radiation, contact resistance,
  electrical displacement current, nonlinear volumetric heat capacity,
  independently certified energy-flux reconstruction or error estimator.
- A converged algebraic fixed point does not establish solution uniqueness,
  continuum accuracy, global nonlinear stability or commercial readiness.

Keep this feature draft until exact-head Ubuntu, Windows, installed-SDK and
Linux address/leak/undefined-behavior sanitizer CI all pass.
