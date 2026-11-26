# Maxwell 3D cell result export

Use `nova/io/vtk_maxwell_3d.hpp` with the driven isotropic H(curl) Maxwell SDK:

```cpp
auto fields = nova::io::reconstruct_maxwell_cells_3d(
    mesh, materials, k0, system, solution.edge_field);
nova::io::write_legacy_vtk_maxwell_3d(output, mesh, fields);
```

Check solver convergence before export. E is affine within each tetrahedron; H is constant, so centroid E, H and time-averaged Poynting flux are exact cell averages. Loss density is quadratic: the mean uses exact cell integration divided by cell volume, not centroid loss. Complex fields are peak phasors with exp(-i omega t).

All output is CELL_DATA (one tuple per element). electric_real/imag are V/m; magnetic_real/imag are A/m; poynting_mean is W/m^2; absorbed_power is W; mean_absorbed_power_density is W/m^3. Separate real/imag fields preserve phase. ParaView can open the legacy ASCII VTK output. This piecewise-constant visualization does not reconstruct the full affine field or imply nodal continuity across material interfaces. No new solver, port, PML, anisotropic or nonlinear support is added.

Regression cases cover permuted tetrahedron affine fields and exact loss, multi-element constant field and zero loss, malformed edge ordering, invalid loss, stream failure and formatting. CI execution required before merge.
