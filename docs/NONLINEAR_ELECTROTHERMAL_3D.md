# Experimental 3D nonlinear stationary electrothermal feedback

Public API: `nova/physics/nonlinear_electrothermal_3d.hpp`,
`solve_nonlinear_electrothermal_3d`.

This builds on the separately merged, one-way electrical/Joule/thermal FEM
workflow. It is a **stationary temperature-dependent electrical and optional thermal
conductivity** Picard iteration. It is not a transient or fully monolithic multiphysics
solver.

## Model and units

For each electrical P1 tetrahedron `e`, the reference conductivity is
a positive isotropic cell scalar or a symmetric-positive-definite cell tensor.
At iteration `n`:

- `Tbar_e^n = (T_i^n + T_j^n + T_k^n + T_l^n)/4` kelvin.
- `sigma_e^n = sigma_ref_e * scale(Tbar_e^n)`, where `scale` is
  a caller-provided, finite, **strictly positive**, dimensionless law.
- Solve `-div(sigma^n grad(V^n)) = s_I`, with electrical boundary data.
- Calculate `E=-grad(V)`, `J=sigma^n E`, and `q_J=E dot J >= 0`
  using the existing one-way flux/heat assembly.
- If a thermal law is supplied, form `k_e^n = k_ref_e * k_scale(Tbar_e^n)`
  using the same cell-averaged iterate, then solve
  `-div(k^n grad(T_solved^n)) = q_base + q_J`.
- Without a thermal law, the base thermal conductivity remains unchanged,
  including its original spatial-scalar and spatial-tensor representations.
- If not converged, set `T^(n+1) = T^n + relaxation*(T_solved^n-T^n)`.

The two physics must use **identical** tetrahedral connectivity and coordinates.
Electrical conductivity must be provided entirely by
`electrical.cell_diffusivity` or `electrical.cell_diffusivity_tensors`.
Spatial callbacks are intentionally rejected because this API samples a
cell-average nodal temperature, not the quadrature-point temperature.
For temperature-dependent `k(T)`, thermal material must also use precisely
one complete cellwise scalar or SPD-tensor representation (no spatial callback).
The scaling laws must return **strictly positive finite** multipliers.
Without `k(T)`, thermal spatial scalar/tensor/cellwise materials remain
available. Supported Dirichlet, Neumann and Robin conditions are unchanged.

The solver checks **unrelaxed** differences, so tiny relaxation cannot
misleadingly declare success:

```text
max_i |T_solved_i - T_iter_i| /
       (abs_tol_K + rel_tol * max(|T_solved_i|, |T_iter_i|)) <= 1
```

Both tolerances must be finite and nonnegative, and at least one positive.
The relaxation is in (0,1], and the iteration cap is positive. Invalid inputs
and failed inner linear solves throw. Reaching the cap returns a result with
`converged=false`, the final fields, the last temperature iterate, the raw
maximum fixed-point defect and its normalized measure. Check `converged`
before using results for scientific interpretation.

## Optional adaptive relaxation and iteration history

Default behavior remains **fixed relaxation** (0.7), preserving existing
solves. When `options.adaptive_relaxation=true`, the solver compares the
unrelaxed, tolerance-normalized Picard fixed-point defect with the preceding
iteration:

- If the current residual exceeds **105%** of the previous residual, halve
  the relaxation, not below `minimum_relaxation`.
- If the residual falls below **80%** of the previous residual, increase
  relaxation by 20%, not above the original `options.relaxation`.
- Otherwise retain the current relaxation.

The damping decision is made **after** the current full electrical/thermal
solve and affects the update to the next iterate. It does not attempt an
extra trial solve, and does **not guarantee monotonic residual reduction**,
convergence, uniqueness, or a thermodynamically admissible solution.
Adaptive minimum relaxation must be finite, positive and no greater than
the initial `relaxation`.

Each completed paired solve appends one
`NonlinearElectrothermal3DIteration` record to `result.history` with:
iteration number, relaxation applied *after* that solve (0 on convergence
or iteration cap), max absolute fixed-point residual in kelvin, normalized
residual, generated Joule watts, discrete weak-form thermal-balance defect
in watts, and minimum/maximum solved nodal temperatures in kelvin.
`result.history.size() == result.iterations` on normal return.
On nonconvergence, the full history is returned with `converged=false`.

Example:

```cpp
nova::physics::NonlinearElectrothermal3DOptions opts;
opts.max_iterations = 60;
opts.relaxation = 0.7;
opts.adaptive_relaxation = true;
opts.minimum_relaxation = 0.1;
auto result = nova::physics::solve_nonlinear_electrothermal_3d(p, opts);
for (const auto& step : result.history) {
    // step.iteration, step.max_normalized_residual,
    // step.relaxation_applied, step.total_joule_power_watts
}
```

## C++ usage

```cpp
#include "nova/physics/nonlinear_electrothermal_3d.hpp"

nova::physics::NonlinearElectrothermal3DProblem p;
p.base = existing_one_way_problem; // electrical + thermal, identical mesh
// The electrical reference material must be a complete per-cell array:
p.base.electrical.diffusivity = {};
p.base.electrical.cell_diffusivity.assign(
    p.base.electrical.mesh.elements.size(), 2.0); // S/m
p.initial_temperature_kelvin = 300.0;
p.conductivity_scale = [](double t_kelvin) {
    return 1.0 + 0.003 * (t_kelvin - 300.0);
};
// Optional: temperature-dependent thermal conductivity in W/(m K).
// With this enabled, base.thermal must have complete per-cell material.
p.base.thermal.diffusivity = {};
p.base.thermal.cell_diffusivity.assign(
    p.base.thermal.mesh.elements.size(), 4.0);
p.thermal_conductivity_scale = [](double t_kelvin) {
    return 1.0 + 0.01 * (t_kelvin - 300.0);
};
const auto result = nova::physics::solve_nonlinear_electrothermal_3d(p);
if (!result.converged) {
    // Inspect result.iterations and result.max_fixed_point_residual_kelvin.
}
```

## Verification

`tests/test_nonlinear_electrothermal_3d.cpp` covers:

- Zero-feedback equivalence to the established one-way plate benchmark
  (200 W for 2 S/m, 10 V across a unit cube, and 4 W/(m K)).
- Positive feedback, conductivity consistency and voltage gauge invariance.
- Constant-voltage SPD tensor support with additive 10 W base heating.
- Thermal `k(T)` regression showing reduced interior temperature for a
  positive thermal-conductivity slope, with unchanged 200 W Joule power.
- Coupled `sigma(T)` and `k(T)` gauge invariance and energy budget.
- Tensor thermal `k(T)` coefficients with additive spatial+cell heat sources.
- Backward-compatible thermal spatial callback when `k(T)` is absent.
- Invalid thermal scaling laws, mixed thermal materials and incomplete arrays.
- Iteration-cap nonconvergence and tiny-relaxation false-convergence guard.
- Bounded adaptive damping, iteration-by-iteration telemetry and installed-SDK
  access to the exported diagnostics.
- Invalid laws, material representations, data lengths, initial conditions,
  and nonlinear solver options.
- `tests/package_consumer/main.cpp` links the public symbol from the
  installed SDK.

Run the project's Ubuntu/Windows build, CTest, installed-SDK consumer and
Linux ASan/LSan/UBSan jobs on the **exact branch head** before marking the PR
ready for review. Tests and CI status are not asserted as passing until
verified in GitHub Actions.

## Current limitations

- Conductivity is held **constant inside each cell**, based on the cell
  average of the P1 nodal temperature. This is a discretization choice;
  finer cellwise material resolution can change the result.
- Electrical and thermal conductivity can depend on **cell-average** nodal
  temperature through independent dimensionless positive scaling laws.
  Other coefficients and constitutive effects are temperature-independent.
- The returned energy budget is an assembled discrete weak-form balance; it
  is not an independently certified continuum heat-flux calculation.
- Fixed-point convergence is **not proof of uniqueness**, continuum accuracy,
  energy-norm convergence or global nonlinear stability. No automatic line
  search, Newton Jacobian, temperature-dependent contact resistance, radiation,
  Seebeck or Peltier coupling.
- No mesh-to-mesh transfer, different physics meshes, transient thermal
  feedback, automatic continuation, or certified error estimate.
