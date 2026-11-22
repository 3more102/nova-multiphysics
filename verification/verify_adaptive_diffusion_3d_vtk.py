#!/usr/bin/env python3
"""End-to-end check of the adaptive 3D demo's legacy VTK outputs.

Uses only the Python standard library and independently validates the
exported mesh topology, solution values, and volume conservation.
"""
import math
import pathlib
import subprocess
import sys
import tempfile


def read_vtk(path: pathlib.Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    if lines[:4] != [
        "# vtk DataFile Version 3.0",
        "NOVA P1 tetrahedral scalar result",
        "ASCII",
        "DATASET UNSTRUCTURED_GRID",
    ]:
        raise AssertionError(f"Unexpected VTK header in {path}")
    points_header = lines[4].split()
    if len(points_header) != 3 or points_header[0] != "POINTS":
        raise AssertionError("Missing POINTS section")
    n_points = int(points_header[1])
    if n_points == 0:
        raise AssertionError("Export contains zero points")
    points = [
        tuple(float(x) for x in lines[5 + i].split())
        for i in range(n_points)
    ]
    if any(len(point) != 3 or not all(map(math.isfinite, point))
           for point in points):
        raise AssertionError("Invalid point coordinates")

    offset = 5 + n_points
    cells_header = lines[offset].split()
    if len(cells_header) != 3 or cells_header[0] != "CELLS":
        raise AssertionError("Missing CELLS section")
    n_cells, packed = map(int, cells_header[1:])
    if not n_cells or packed != 5 * n_cells:
        raise AssertionError("Invalid tetrahedral connectivity metadata")
    cells = []
    for i in range(n_cells):
        indices = [int(v) for v in lines[offset + 1 + i].split()]
        if len(indices) != 5 or indices[0] != 4:
            raise AssertionError("Invalid tetrahedral cell record")
        cell = indices[1:]
        if any(index < 0 or index >= n_points for index in cell):
            raise AssertionError("Out-of-range VTK vertex index")
        if len(set(cell)) != 4:
            raise AssertionError("Repeated tetrahedron vertex")
        cells.append(cell)

    offset += 1 + n_cells
    if lines[offset] != f"CELL_TYPES {n_cells}":
        raise AssertionError("Wrong cell-type count")
    if lines[offset + 1: offset + 1 + n_cells] != ["10"] * n_cells:
        raise AssertionError("Non-tetrahedral VTK cell type")

    offset += 1 + n_cells
    if lines[offset:offset + 3] != [
        f"POINT_DATA {n_points}",
        "SCALARS temperature double 1",
        "LOOKUP_TABLE default",
    ]:
        raise AssertionError("Missing nodal scalar data")
    values = [float(value) for value in lines[offset + 3:]]
    if len(values) != n_points or any(not math.isfinite(v) for v in values):
        raise AssertionError("Invalid VTK nodal field")

    volume = 0.0
    for cell in cells:
        a, b, c, d = (points[i] for i in cell)
        u = [b[i] - a[i] for i in range(3)]
        v = [c[i] - a[i] for i in range(3)]
        w = [d[i] - a[i] for i in range(3)]
        det = (
            u[0] * (v[1] * w[2] - v[2] * w[1])
            - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0])
        )
        volume += abs(det) / 6.0
    if abs(volume - 1.0) > 1.0e-10:
        raise AssertionError(f"Mesh volume changed: {volume}")

    for point, value in zip(points, values):
        if point[0] in (0.0, 1.0) and abs(value) > 1.0e-8:
            raise AssertionError("Dirichlet boundary not preserved in VTK")
        if value < -1.0e-8 or value > 0.26:
            raise AssertionError("Unphysical manufactured wall temperature")
    return n_points, n_cells


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: verify_adaptive_diffusion_3d_vtk.py executable")
    with tempfile.TemporaryDirectory(prefix="nova_adaptive3d_") as temp:
        directory = pathlib.Path(temp)
        coarse = directory / "coarse.vtk"
        refined = directory / "refined.vtk"
        process = subprocess.run(
            [sys.argv[1], "--vtk", str(coarse), str(refined)],
            capture_output=True, text=True, check=True, timeout=120,
        )
        if "NOVA adaptive 3D diffusion demonstration" not in process.stdout:
            raise AssertionError("Missing expected example diagnostics")
        n0, k0 = read_vtk(coarse)
        n1, k1 = read_vtk(refined)
        if n1 <= n0 or k1 <= k0:
            raise AssertionError("Adaptive mesh did not increase resolution")
        print(f"Adaptive 3D VTK verified: nodes {n0}->{n1}, cells {k0}->{k1}")


if __name__ == "__main__":
    main()
