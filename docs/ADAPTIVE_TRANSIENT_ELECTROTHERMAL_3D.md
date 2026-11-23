# Experimental adaptive 3D transient electrothermal integration

Public C++20 API:
`nova/physics/adaptive_transient_electrothermal_3d.hpp`,
`solve_adaptive_transient_electrothermal_3d`.

This builds on the merged [implicit Backward Euler coupled solver](TRANSIENT_ELECTROTHERMAL_3D.md)
and its quasi-static DC electrical field, P1 Joule source, thermal capacity,
temperature-dependent sigma(T) and optional k(T). Spatial mesh remains fixed.

## Temporal error control

At each physical time t_n, the controller constructs **three independent**
nonlinear solves, starting from the same last accepted thermal state:

1. One full Backward Euler step, length h.
2. A first half-step, h/2, starting at t_n.
3. A second half-step, h/2, starting from **accepted first-half candidate**
   and ending at t_n+h.

The controller estimates the local time-discretization error using the
**unextrapolated, full-vs-two-half P1 nodal difference**:

```
E = max_i abs(T_half2[i] - T_full[i]) /
           (abs_error_K + rel_error * max(abs(T_half2[i]),abs(T_full[i])))
```

This comparison has leading O(h^2) behavior for smooth Backward Euler
solutions (first-order global accuracy). It is an *indicator*, not a
certified local or continuum error bound.

The candidate is accepted only if **all three** segregated nonlinear
solves converge and E <= 1. On acceptance, the final state is the **two
half-step** result, never the less-accurate full-step solution.
On rejection, none of the three candidate temperatures is committed.

For successful nonlinear comparisons, the new trial step is chosen by
`safety_factor / sqrt(E)` (bounded by a shrink factor and growth factor).
When E=0 the growth factor is used. Rejections always shrink (at most 90%
of the current h). If a nonlinear fixed-point trial fails, the controller
shrinks by one half (or the configured shrink floor). A requested final
interval shorter than `minimum_time_step` may be tried to reach
`end_time`, but a rejected interval at/below minimum stops safely.

If `maximum_attempts`, `maximum_accepted_steps`, `minimum_time_step`
or representable time resolution prevents completion, the function returns
`completed=false` with an explicit `stop_reason`. The returned
`temperature` and `final_time` always refer to the **last accepted**
state; no failed trial is committed. Exceptions from invalid inputs,
invalid constitutive laws, and failed inner linear solvers still propagate.

## Thermal and electrical energy diagnostics

Every comparison recalculates a DC field and cellwise Joule heat using
the thermal iterate of that candidate. Its energy units are:

- `joule_power_watts`: W at the accepted second-half-step endpoint.
- `heat_budget.source_energy`: J added over both accepted half-steps.
- `heat_budget.heat_content`: J evaluated at accepted endpoint.
- Neumann, Robin and Dirichlet energies: **sums of both half-steps**, J.
- `heat_budget.balance_residual`: sum of both half-step algebraic
  weak-form energy-balance defects in J.

The rejected **full** step's heat budget is never added to accepted
energy accounting. `attempts` logs start, size, end, normalized error,
nonlinear success and acceptance. `accepted_steps` records physical dt,
error, nonlinear iteration cost for three solves, Joule power and energy.
`rejected_steps` counts attempts not accepted. The nonlinear iteration
count includes the discarded full-step comparison (cost diagnostic).

## Example

```cpp
#include "nova/physics/adaptive_transient_electrothermal_3d.hpp"

nova::physics::AdaptiveTransientElectrothermal3DOptions options;
options.coupling.thermal_time.start_time = 0.0;
options.coupling.thermal_time.scheme =
    nova::fem::TransientHeat3DScheme::BackwardEuler;
options.end_time = 1.0;
options.initial_time_step = 0.1;
options.minimum_time_step = 1.0e-6;
options.maximum_time_step = 0.25;
options.absolute_error_kelvin = 1.0e-4;
options.relative_error = 1.0e-6;
const auto r = nova::physics::solve_adaptive_transient_electrothermal_3d(
    transient_electrothermal_problem, options);
if (!r.completed) {
    // r.stop_reason, r.final_time, r.rejected_steps and r.attempts
}
```

## Numerical verification cases

`tests/test_transient_electrothermal_3d.cpp` exercises:

- Manufactured unit-cube DC field: 10 V across 1 m,
  sigma=2 S/m, C=2 J/(m^3 K), insulated heat -> 200 W
  Joule power and exact T(t)=300+100t; 28 J over 0.14 s.
- Nonlinear uniform feedback sigma(T)=2*[1+0.01*(T-300)]
  yields the independent ODE **d(T-300)/dt=100+(T-300)**.
  Its exact solution is T(t)=300+100*(exp(t)-1). Adaptive
  integration rejects an overly large first step and improves on
  single-step Backward Euler at t=0.3 s, including k(T) scaling.
- Nonlinear Picard failure at max iterations=1 causes successive
  rejected steps, reaches minimum dt, and returns the unchanged
  300 K initial state without silently advancing time.
- Independent attempt/accepted-step limits, short final subminimum
  interval, invalid step controls, input immutability, energy accounting.
- Installed SDK consumer compiles, links and uses public step history.

## Limitations

- No spatial mesh refinement, temperature-dependent capacity, stiffness
  Jacobian, transient electrical displacement or multi-rate coupling.
- No embedded higher-order time method: full-vs-half Backward Euler
  controls only the **temporal** discrepancy of nonlinear FEM solutions.
- Three distinct electrothermal nonlinear solves are required for
  every comparison, so high-accuracy simulation may be computationally
  expensive.
- Underlying per-time-step nonlinear tolerance must be sufficiently
  tighter than requested temporal error tolerance; otherwise local error
  measurement may be polluted by algebraic solve error.
- Endpoint time-dependent boundary conditions are evaluated through the
  existing transient solver; rejected trial callbacks can be re-evaluated
  during retries. Callbacks should be pure and deterministic.
- Strong nonlinear instability may cause minimum-step termination;
  automatic controller success is not a proof of nonlinear uniqueness
  or physical fidelity.
- Experimental numerical code; not commercially validated.

Run exact-head Ubuntu, Windows, installed SDK and Linux
address/leak/undefined-behavior sanitizer CI before declaring this
feature validated or merging.
