# Adaptive transient 2D heat (experimental)

`nova/fem/adaptive_heat_2d.hpp` provides fixed-mesh temporal step doubling
for the linear P1/P2 H1 heat equation. Each attempted interval compares one
full step with two half steps from the same accepted DOF vector. The
two-half-step solution is accepted when the normalized RMS error is at most
one; rejected attempts do not change the accepted state or energy history.

## Prepared operators and cache ownership

At the start of an adaptive study NOVA assembles the spatially varying,
time-independent capacity and conductivity operators once. It also prepares
the H1 boundary-DOF map and capacity weights. Separate full-step and
half-step cache slots retain their respective time-step-dependent base
operator (M + \theta \Delta t K).

The source, Dirichlet trace, Neumann flux, Robin transfer coefficient and
ambient temperature remain time-dependent. They are evaluated at the
required physical time levels for each trial. The current Robin matrix is
added to the cached base operator before Dirichlet elimination. Consequently,
the complete constrained matrix can change even when the time step does not.

Each slot owns a `CachedConjugateGradient`. IC(0) is reused only when the
complete constrained CSR matrix is exactly identical; a changed time step,
Robin coefficient or sparsity/coefficient value causes a rebuild. The result
reports:

- `full_step_base_operator_builds` and `half_step_base_operator_builds`.
- `full_step_ic0_factor_builds` and `half_step_ic0_factor_builds`.

The base-operator counters measure distinct time-step/scheme-dependent
(M + \theta \Delta t K) preparations, not the time-dependent Robin
assembly. IC(0) counters are zero when a different preconditioner is selected.

## Verification

`tests/test_adaptive_heat_2d.cpp` covers P1/P2 constant-state preservation,
time-dependent source and Dirichlet data, adaptive rejection and temporal
accuracy, energy-budget closure, fixed-step restart, one-time material
operator preparation, repeated-`dt` base-operator reuse, time-varying Robin
matrix invalidation, and exact IC(0) reuse.

## Scope

This is fixed-mesh temporal adaptivity, not spatial mesh refinement. Material
capacity and conductivity remain spatially varying but time-independent.
The method does not provide a global temporal error bound or engineering
certification; NOVA remains an experimental numerical platform.
