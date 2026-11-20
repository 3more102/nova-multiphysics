# Validated 3D thermal restart checkpoints (experimental)

NOVA's fixed-step P1 tetrahedral heat solver already supports restarting from an
absolute nodal-temperature vector at any finite physical start time. This
module adds **versioned text serialization** for that state without changing
the solver or the physics formulation.

## Public API

`#include "nova/io/thermal_checkpoint_3d.hpp"`

```cpp
// The caller owns the stream, mesh, problem and time-stepping settings.
std::ofstream checkpoint_file("thermal.restart");
nova::io::write_transient_heat_checkpoint_3d(
    checkpoint_file, problem.diffusion.mesh,
    first_result.final_time, first_result.temperature);
checkpoint_file.close();
if (!checkpoint_file) throw std::runtime_error("Cannot finish restart file");

// In a new process, reconstruct the same model first.
std::ifstream source("thermal.restart");
const auto saved = nova::io::read_transient_heat_checkpoint_3d(
    source, problem.diffusion.mesh);
problem.initial_temperature = saved.temperature;
options.start_time = saved.time;
const auto second_result = nova::fem::solve_transient_heat_3d(problem, options);
```

Include `<fstream>` and `<stdexcept>` in this example's translation unit.

## Format and safety guarantees

Version 1 is UTF-8-compatible, whitespace-separated ASCII:

```text
NOVA_THERMAL_CHECKPOINT_3D 1
time <double>
nodes <count>
<x> <y> <z>    # repeated for every node
elements <count>
<n0> <n1> <n2> <n3>    # repeated for every tetrahedron
temperature <count>
<T>           # one absolute temperature per node
end
```

Numbers are written with `std::numeric_limits<double>::max_digits10` using
the classic decimal locale, allowing exact finite-double roundtrips. The
reader compares node coordinates and tetrahedron vertex ordering **exactly**
to the already-validated expected mesh; a similarly sized mesh with different
geometry or connectivity is not accepted. It rejects unsupported versions,
truncated/malformed input, non-finite time or temperature, inconsistent counts
and unexpected trailing tokens. The writer validates inputs before emitting
bytes and does not change the caller's stream precision or flags.

The file is **not** a self-contained simulation model: it does not store
conductivity/capacity, sources, boundary conditions, integration scheme, time
step, solver tolerances or mesh physical tags. The caller must reconstruct
the same physical model and numerical scheme when continuing a study. Passing
mesh validation alone does not establish model equivalence.

Writing through a stream does not guarantee atomic file replacement if the
disk becomes full or an I/O error occurs. Applications should write to a new
temporary file, check close/flush success and atomically rename it when
supported by the filesystem. Text output is currently staged in memory to
avoid partial files from **input validation** failures.

The module intentionally does not add checkpoint CLI flags or support mesh
adaptation, project serialization, MPI or cross-platform binary formats.

## Verification

`tests/test_thermal_checkpoint_3d.cpp` checks exact-double roundtrip,
stream formatting, malformed/truncated/trailing input and version rejection,
non-finite value rejection without partial serialization, mesh mismatch
detection (including same-size changes), and full-vs-resumed 3D heat
integration for Backward Euler and Crank–Nicolson.

This validates discrete restart equivalence, not engineering accuracy or
experimental calibration.
