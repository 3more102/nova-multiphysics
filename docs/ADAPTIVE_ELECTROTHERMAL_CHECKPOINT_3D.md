# Versioned adaptive electrothermal restart (experimental)

The single-stream C++20 checkpoint API in
`nova/io/adaptive_electrothermal_checkpoint_3d.hpp` stores an **accepted**
3D electrothermal step boundary as one text record. This replaces the
earlier, application-managed pair of a thermal checkpoint and a separately
stored proposed controller timestep.

The record includes:

- the accepted physical time and full nodal P1 temperature vector;
- `next_time_step`, the adaptive controller's proposed next trial size;
- exact mesh node coordinates, tetrahedron indices and node ordering;
- an explicit, caller-assigned `model_id` compatibility token;
- the remaining final time, adaptive error tolerances/controller factors,
  attempt/accepted-step caps, nonlinear iteration/tolerances and
  electrical and thermal linear-solver numerical settings.

The embedded thermal state is the existing
[`NOVA_THERMAL_CHECKPOINT_3D`](THERMAL_CHECKPOINT_3D.md)
version 1 format. Both layers use lossless finite-double serialization,
mesh checking, strict field names, version rejection and trailing-data
rejection. A failure to validate the input leaves the output stream
untouched; a failed physical stream write can still leave a partial file.

## Public API example

```cpp
#include "nova/io/adaptive_electrothermal_checkpoint_3d.hpp"
#include "nova/physics/adaptive_transient_electrothermal_3d.hpp"

// Recreate the same physical model in both processes.
constexpr auto model_id = "my-electrothermal-study-v1";
nova::physics::AdaptiveTransientElectrothermal3DOptions options;
options.end_time = 1.0;

// Called only after an accepted adaptive step:
options.on_accepted_step = [&](const auto& step, const auto& field) {
    nova::io::write_adaptive_electrothermal_checkpoint_3d_file_atomic(
        "checkpoint.dat", problem.electrical.mesh, model_id, options,
        {step.time, field, step.next_time_step});
};

// Load before changing the starting time or proposed next dt:
const auto state = nova::io::read_adaptive_electrothermal_checkpoint_3d_file(
    "checkpoint.dat", problem.electrical.mesh, model_id, options);
problem.thermal.initial_temperature = state.temperature;
options.coupling.thermal_time.start_time = state.time;
options.initial_time_step = state.next_time_step;
const auto resumed =
    nova::physics::solve_adaptive_transient_electrothermal_3d(
        problem, options);
```

The example assumes `problem` is a configured
`TransientElectrothermal3DProblem`. The original stream APIs remain
available. The accepted-step callback stays registered after restart.

**Load before changing** `start_time` or `initial_time_step`. These two
values are dynamic and intentionally not compared against the record.
All other stored numerical settings must match the read-time options
exactly. The checkpoint's time must fall within the configured physical
window, and the proposed timestep must stay within the configured bounds.

## Format version 1

```text
NOVA_ADAPTIVE_ELECTROTHERMAL_CHECKPOINT_3D
version 1
model_id <opaque-identifier>
end_time <double>
min_dt <double>
max_dt <double>
absolute_error_kelvin <double>
relative_error <double>
safety_factor <double>
minimum_shrink_factor <double>
maximum_growth_factor <double>
maximum_attempts <unsigned>
maximum_accepted_steps <unsigned>
coupling_max_iterations <unsigned>
coupling_relaxation <double>
coupling_abs_kelvin <double>
coupling_rel_tolerance <double>
thermal_max_iterations <unsigned>
thermal_relative_tolerance <double>
thermal_absolute_tolerance <double>
thermal_preconditioner <enum-integer>
electrical_max_iterations <unsigned>
electrical_relative_tolerance <double>
electrical_absolute_tolerance <double>
electrical_preconditioner <enum-integer>
next_dt <double>
thermal_state
NOVA_THERMAL_CHECKPOINT_3D 1
... validated mesh / time / P1 field ...
end
```

The model ID must be 1–128 ASCII characters from `A–Z a–z 0–9 _ - . :`.
It is **not** a content-derived hash, signature or attestation. Use a new
identifier when the material models, heat sources, boundary conditions,
mesh tags, electrical or thermal problem definitions change.

## Boundaries and scientific interpretation

This is **not a fully self-contained physical simulation archive**.
The formulation, constitutive functions, source/BC callbacks,
geometric region semantics, current electrical field and prior
accepted/rejected attempts are not serialized. Mesh equality, ID match and
numerical-settings equality are necessary but not sufficient evidence
of physical model equivalence. Reconstruct the same model and use
deterministic callbacks to obtain reproducible step-boundary restart.

Attempt and accepted-step limits apply to **each API invocation**, not
to all previous invocations combined. Separate application bookkeeping
is required for cumulative compute-budget auditing and history/energy
accounting. No MPI state, adaptive mesh genealogy, time-dependent
external solver state or implicit state in callbacks is captured.

The stream writer stages serialization in memory only. The new
`write_adaptive_electrothermal_checkpoint_3d_file_atomic` path API
serializes the same version-1 record, reserves a unique sibling staging
directory, writes and closes a file, and renames it over the destination.
Failed validation, write, or rename operations leave an existing destination
unchanged, and normal failures clean up staging files. The path-level
`read_adaptive_electrothermal_checkpoint_3d_file` uses the same strict parser.

**Filesystem limits:** Replacement atomicity depends on host filesystem
rename semantics. The implementation does not fsync file data or parent
directory metadata, so it is **not a power-loss durability guarantee**.
Abrupt process death may leave stale `.nova-stage-*` directories beside
the destination; remove only when no writer is active. Concurrent writers
targeting the same destination are not coordinated. The implementation
still does not persist cumulative attempt counters or accepted energy
history, and does not reconstruct material/source callbacks. Treat any
uncertain filesystem failure as a failed publish.

## Verification

The dedicated tests validate exact-double roundtrip, expected-model and
controller incompatibility, altered mesh ordering/coordinates, version
corruption, trailing/truncated records, invalid writes without partial
stream output and a nonlinear electrothermal full-vs-restarted integration.
The installed C++20 SDK consumer exercises both stream and path-level symbols.

These are discrete reproducibility tests, not continuum accuracy
certification or a demonstration of COMSOL/ANSYS equivalence.
