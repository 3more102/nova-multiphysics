# ParaView time-series PVD collections

`nova::io::write_vtk_pvd_collection` writes a VTK XML Collection (`.pvd`) file that
groups existing `.vtk` or `.vtu` snapshots into a time series for ParaView.
It does **not** create the snapshot files. The existing
`write_legacy_vtk_scalar_3d` writes the P1 tetrahedral temperature snapshots.

## Example with the 3D transient heat solver

```cpp
#include "nova/fem/transient_heat_3d.hpp"
#include "nova/io/vtk_legacy_3d.hpp"
#include "nova/io/vtk_pvd.hpp"
#include <fstream>
#include <string>
#include <vector>

// Here, problem and options are already configured.
const auto result = nova::fem::solve_transient_heat_3d(problem, options);
std::vector<nova::io::VtkTimeSeriesEntry> frames;
for (std::size_t i = 0; i < result.linear_solvers.size(); ++i) {
    const auto filename = "temperature_" + std::to_string(i + 1) + ".vtk";
    std::ofstream vtk(filename);
    nova::io::write_legacy_vtk_scalar_3d(
        vtk, problem.diffusion.mesh, result.linear_solvers[i].solution);
    frames.push_back({result.diagnostics[i].time, filename});
}
std::ofstream collection("temperature.pvd");
nova::io::write_vtk_pvd_collection(collection, frames);
```

Open `temperature.pvd` in ParaView and use its time controls. The example
starts at the first completed step; include a separately exported initial
snapshot if the initial condition is also required in the animation.

## Validation and portability

- Times are finite and strictly increasing; negative start times are allowed.
- Each dataset filename is unique, relative, ends in `.vtk` or `.vtu`,
  and contains only portable printable ASCII path characters.
- Absolute paths, `.`/`..` segments, empty segments, backslashes, and
  unsupported path characters are rejected to avoid ambiguous references.
- XML attribute characters such as `&` and apostrophes are escaped.
- Time is printed with `max_digits10` precision in the C locale.
- All inputs are checked before writing. A failed output stream raises an exception.

Paths are interpreted relative to the directory of the `.pvd` file. The
caller remains responsible for creating snapshot files at those paths, for
ensuring their content matches the listed physical time, and for creating
any needed directories. This is a manifest feature, not a VTK data reader
or an independent physical validation of the transient solver.
