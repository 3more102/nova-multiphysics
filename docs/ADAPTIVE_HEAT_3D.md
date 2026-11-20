# Adaptive transient 3D heat (experimental)

`nova/fem/adaptive_heat_3d.hpp` adds **fixed-mesh temporal adaptivity** to
NOVA's consistent-mass P1 tetrahedral transient solver. Conductivity and
volumetric capacity may be static or prescribed time-dependent scalar/tensor
fields. Time-dependent volumetric sources, Dirichlet nodal/face values,
Neumann face fluxes and Robin transfer/ambient data are sampled consistently
by both full-step and half-step trials.

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
// Optional accepted-state stream: index 0 is the projected initial field.
options.on_step = [](std::size_t index, double time,
                     const std::vector<double>& temperature) {
    // Persist or visualize this field while the callback is active.
};
const auto result = nova::fem::solve_adaptive_heat_3d(problem, options);
```

`TransientHeat3DProblem::time_dependent_source(x,y,z,t)` adds to the existing
static spatial and cell sources. Backward Euler evaluates the source at each
substep's new time; Crank–Nicolson averages the source vectors at both endpoints.
The full-step estimate and the two half-step estimate therefore sample the
same physical source history at their respective time grids.

Time-dependent scalar conductivity, SPD tensor conductivity and volumetric
capacity callbacks are evaluated at the start and end of each trial substep.
Backward Euler uses the old capacity matrix on the old state and the new
material operators implicitly; Crank–Nicolson uses both endpoint matrices.
When capacity varies, storage is conservative, `d(C*T)/dt`, and the heat
content at each endpoint is evaluated with that endpoint's capacity matrix.
The full-step and two-half-step estimates therefore use identical material
time-level semantics.

Time-dependent Dirichlet constraints and faces are evaluated at the study's
physical start time and at each substep endpoint. Time-dependent Neumann face
fluxes are evaluated at the new endpoint for Backward Euler and at both
endpoints for Crank–Nicolson, using trapezoidal load integration. Time-varying
Robin transfer and ambient callbacks are also evaluated at each required
endpoint. Ambient-only changes affect the right-hand side; changes to the
transfer coefficient alter the Robin matrix and trigger exact constrained-CSR
cache comparison and IC(0) rebuild when necessary.


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
accepted-only observer ordering, initial Dirichlet projection, observer
exception propagation, repeated-dt operator/factor reuse and cache rebuilds
when the controller changes the time step. `tests/test_transient_heat_3d.cpp`
also verifies time-dependent Robin endpoint matrices, dynamic scalar/tensor
conductivity, dynamic capacity, conservative energy closure and BE/CN temporal
convergence.

## Assembly and factor reuse

The adaptive study prepares validated mesh topology, static material
operators, static loads and Dirichlet data once. Dynamic conductivity and
capacity matrices are integrated at the material time levels of each full or
half-step trial. Time-dependent source and Neumann loads are evaluated at the
temporal quadrature points; dynamic Robin matrices and ambient loads are
evaluated at the same physical endpoint times. Separate full-step and
half-step caches compare the previous/next mass and conductivity matrices,
Robin matrices, and `dt` exactly before reusing constrained CSR operators.
A changed material matrix rebuilds the affected operator; IC(0) is rebuilt
only if the complete constrained CSR system changes. The result reports
`full_step_operator_builds`, `half_step_operator_builds`,
`full_step_ic0_factor_builds` and `half_step_ic0_factor_builds` so applications
and regression tests can observe this behavior. IC(0) build counters remain
zero when another preconditioner is selected.

The cache retains only the most recent time step in each slot, not a map of all
previous step sizes. A later return to an older `dt` after an intervening change
therefore rebuilds that slot's operators and, when required, its IC(0) factor.
This is a performance optimization, not a change to the temporal error
estimator or acceptance rule. NOVA remains experimental and is not
engineering-design certification.


## Spatial adaptation building blocks

The 3D mesh layer now provides deterministic Dörfler bulk marking and
one-level conforming local tetrahedral refinement. Given nonnegative
element indicators, mark the smallest descending-indicator set whose squared
sum reaches the requested fraction of the global squared sum. The marking
routine rejects non-finite or negative indicators and requires a fraction in
(0, 1].

Each marked tetrahedron contributes its six edges to a global split set.
All incident tetrahedra are bisected along those same marked edges in a
deterministic order, so both sides of every shared face construct matching
subdivisions without hanging faces. The returned refinement map records each
child's parent, each new node's coarse edge, and each refined exterior
boundary face's parent. Generic P1 and cellwise-constant transfer routines
consume these maps.

These APIs are building blocks, not an automatic spatially adaptive heat
solver. The current workflow does not yet perform repeated estimate-mark-
refine-resolve cycles, multilevel grading, shape-quality optimization,
physical-group remapping, or automatic reconstruction of material and
boundary callbacks. Applications must transfer the P1 state, map cellwise
data, preserve their physical boundary metadata, and inspect the refined
mesh quality before solving again.
