# Adaptive transient 3D heat (experimental)

`nova/fem/adaptive_heat_3d.hpp` adds **fixed-mesh temporal adaptivity** to
NOVA's existing consistent-mass P1 tetrahedral transient solver. It does not
change the supported materials or boundary conditions: coefficients and
Dirichlet/Neumann/Robin data remain **time-independent**.

## API

```cpp
#include "nova/fem/adaptive_heat_3d.hpp"

nova::fem::AdaptiveHeat3DOptions options;
options.scheme = nova::fem::TransientHeat3DScheme::CrankNicolson;
options.start_time = 0.0;
options.end_time = 1.0;
options.initial_step = 0.1;
options.min_step = 1.0e-5;
options.max_step = 0.2;
options.absolute_tolerance = 1.0e-5;
options.relative_tolerance = 1.0e-4;
const auto result = nova::fem::solve_adaptive_heat_3d(problem, options);

// Optional accepted-state stream: index 0 is the projected initial field.
options.on_step = [](std::size_t index, double time,
                     const std::vector<double>& temperature) {
    // Persist or visualize this field while the callback is active.
};
```

For each attempted interval `dt`, the integrator runs one full step and
two half steps **from the same previously accepted nodal state**. The order
is `p=1` for Backward Euler, `p=2` for Crank–Nicolson. Its RMS nodal
indicator is

```text
e = sqrt(mean_i( [ (T_2half[i] - T_full[i]) / (2^p-1) /
    (atol + rtol * max(abs(T_old[i]), abs(T_2half[i]))) ]^2 ))
```

If `e <= 1`, the two-half-step state is accepted. Otherwise the trial is
rejected without modifying the accepted nodal field or energy history. The
next interval uses `0.9 * e^(-1/(p+1))`, bounded to control growth or
shrinkage. A smaller final interval is allowed to land on `end_time`.

The result includes accepted/rejected counts, a record of all attempted
intervals and their errors, a final nodal field, and **two linear solver
diagnostics for each accepted step** (with full field copies suppressed).
Each accepted energy diagnostic is the sum of both half-step contributions;
`heat_content` is taken from the **second** half step. Thus rejected trial
energy is never accidentally included. Restart by passing the returned
temperature to `TransientHeat3DProblem::initial_temperature` and using
`result.final_time` as the next `start_time`.

`AdaptiveHeat3DOptions::on_step` emits the Dirichlet-projected initial field
once at index 0, then only the accepted two-half-step macro fields at indices
1..N. The coarse full-step estimate, intermediate half-step field and
rejected trials are never emitted. The field reference is valid only during
the callback. Observer exceptions propagate and abort the integration.

## Verification and limits

`tests/test_adaptive_heat_3d.cpp` covers exact insulated uniform heating
and heat budget closure for both schemes, temporal tolerance sensitivity
against the analytical semidiscrete one-free-DOF decay solution, mixed
boundary restart, invalid control rejection, the maximum-attempt gate,
accepted-only observer ordering, initial Dirichlet projection and observer
exception propagation.

**Performance limitation:** each full-step and two-half-step trial invokes
the existing solver, so the mesh, FE operators, and linear systems are
reassembled for each trial. Caching the assembled operators across trial
steps and efficiently reusing factors is a future optimization. This is a
research-grade implementation, not engineering-design certification.
