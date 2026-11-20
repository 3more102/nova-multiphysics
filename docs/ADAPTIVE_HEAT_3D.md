# Experimental Adaptive 3D Heat Time Integration (mainline)

`solve_adaptive_heat_3d` in `nova/fem/adaptive_heat_3d.hpp` adds
fixed-mesh adaptive time stepping to the **current** P1 tetrahedral heat solver.
It is a correctness-first forward port of the concept in older PR #20,
whose prior transient model API is incompatible with current `main`.

For Backward Euler (order 1) or Crank–Nicolson (order 2), each attempt
compares one step of duration `dt` against two consecutive half steps.
The nodal maximum norm of their difference, divided by `2^p-1`, is
scaled by `absolute_tolerance + relative_tolerance * max(|old|,|full|,|fine|)`.
A normalized error of at most 1 accepts the two-half-step state; otherwise
the trial is rejected and the proposed step is reduced. The final step is
clipped to reach `end_time` and may be shorter than `min_step`; a rejected
step at or below the minimum fails explicitly.

Accepted attempts provide an aggregate discrete heat budget: the endpoint
heat content and summed source, Neumann, Robin, Dirichlet-reaction energies,
and balance residuals from both half steps. Rejected attempts never change
the accepted state or published energy diagnostics. `on_step` publishes
the Dirichlet-projected initial field (index 0) and accepted macro-step fields
only. Exceptions propagate. The returned exact nodal field can be supplied
as `initial_temperature` for restart with `start_time=final_time`.

## Example

```cpp
nova::fem::AdaptiveHeat3DOptions opts;
opts.scheme = nova::fem::TransientHeat3DScheme::CrankNicolson;
opts.end_time = 5.0;
opts.initial_step = 0.1;
opts.absolute_tolerance = 1e-6;
opts.relative_tolerance = 1e-4;
auto result = nova::fem::solve_adaptive_heat_3d(problem, opts);
```

## Limitations

This first forward port delegates each full and half step to the existing
validated fixed-step solver; **the mass/stiffness matrices are currently
reassembled for each substep**, and no IC(0) factors are shared across
trials. It is therefore not a performance-optimized adaptive implementation.
Material coefficients, heat sources and boundary values must be
time-independent, as required by the underlying solver. The norm controls
a local step-doubling defect, not a verified global error bound. This is
not a substitute for mesh adaptivity, independent solver verification,
or a commercial-grade thermal certification workflow.

Reproduce: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel && ctest --test-dir build --output-on-failure`.
