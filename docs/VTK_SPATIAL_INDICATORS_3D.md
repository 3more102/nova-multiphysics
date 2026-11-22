# Visualizing tetrahedral spatial error indicators (experimental)

The legacy ASCII VTK writer now supports simultaneous **point-associated**
P1 solutions and **cell-associated** scalar diagnostics on a validated
tetrahedral mesh. Open the resulting `.vtk` file in ParaView and select the
`error_indicator` cell array to inspect where an adaptive FEM algorithm
would place refinement effort.

```cpp
#include "nova/io/vtk_legacy_3d.hpp"

#include <fstream>
#include <vector>

// Assume `mesh`, `nodal_solution` and `element_indicators` are
// independently obtained from the FEM solver and residual estimator.
// element_indicators[k] corresponds exactly to mesh.elements[k].
std::ofstream output("spatial_residuals.vtk");
nova::io::write_legacy_vtk_point_and_cell_scalar_3d(
    output, mesh, nodal_solution, element_indicators,
    "solution", "error_indicator");
```

The writer exports VTK legacy `UNSTRUCTURED_GRID` with linear tetrahedra:
- `POINT_DATA`: one scalar value per mesh vertex.
- `CELL_DATA`: one scalar value per tetrahedron, preserving element ordering.

Each scalar array name must be an ASCII identifier; non-finite or
mismatched-sized data are rejected before any output is written. The writer
uses classic-locale, lossless double-precision text and restores the
caller's output formatting afterward. It detects stream errors but, as with
the existing VTK writer, cannot roll back an I/O stream after a device fails
mid-write.

This export neither computes an error estimate nor proves its reliability.
Residual indicators are local **heuristics**, not certified error bounds;
compare them against manufactured-solution energy errors and mesh-convergence
studies before interpreting a visualization as quantitative accuracy.
