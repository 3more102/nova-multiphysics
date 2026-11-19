#!/usr/bin/env python3
"""End-to-end check that 3D VTK time series matches the monolithic thermal solve.

Only the Python standard library is required; runs inside a temporary folder.
"""
from __future__ import annotations

import math
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def run(executable: Path, *args: str, cwd: Path) -> str:
    result = subprocess.run(
        [str(executable), *args], cwd=cwd, text=True,
        capture_output=True, check=True,
    )
    return result.stdout


def read_values(path: Path) -> list[float]:
    lines = path.read_text(encoding="ascii").splitlines()
    if "DATASET UNSTRUCTURED_GRID" not in lines:
        raise AssertionError(f"{path.name}: not an unstructured VTK grid")
    point_line = next(s for s in lines if s.startswith("POINTS "))
    count = int(point_line.split()[1])
    if f"POINT_DATA {count}" not in lines:
        raise AssertionError(f"{path.name}: incorrect point-data length")
    index = lines.index("LOOKUP_TABLE default")
    values = [float(s) for s in lines[index + 1:]]
    if len(values) != count or not all(map(math.isfinite, values)):
        raise AssertionError(f"{path.name}: wrong length or non-finite data")
    return values


def verify(executable: Path) -> None:
    steps = 4
    with tempfile.TemporaryDirectory(prefix="nova-pvd-check-") as temp:
        folder = Path(temp)
        run(executable, "3", str(steps), "--vtk-series", "case", cwd=folder)
        rows = run(executable, "3", str(steps), "--csv", cwd=folder).splitlines()
        if not rows or rows[0] != "x,y,z,temperature":
            raise AssertionError("CSV baseline missing required header")
        baseline = [tuple(map(float, row.split(","))) for row in rows[1:]]
        collection = folder / "case.pvd"
        root = ET.parse(collection).getroot()
        if root.tag != "VTKFile" or root.attrib.get("type") != "Collection":
            raise AssertionError("PVD is not a VTK Collection")
        datasets = root.findall("./Collection/DataSet")
        if len(datasets) != steps + 1:
            raise AssertionError("PVD does not contain initial and all step states")

        snapshots = []
        for step, dataset in enumerate(datasets):
            expected_name = f"case_{step:06d}.vtk"
            if dataset.attrib.get("file") != expected_name:
                raise AssertionError(f"Unexpected frame filename at {step}")
            expected_time = 0.02 * step / steps
            actual_time = float(dataset.attrib["timestep"])
            if not math.isclose(actual_time, expected_time, abs_tol=1e-14):
                raise AssertionError(f"Incorrect physical time at {step}")
            snapshots.append(read_values(folder / expected_name))

        if len(snapshots[0]) != len(baseline):
            raise AssertionError("Mismatch between VTK and CSV nodal counts")
        for index, (x, y, z, final_value) in enumerate(baseline):
            initial = math.sin(math.pi * x) * math.sin(math.pi * y) * math.sin(math.pi * z)
            if not math.isclose(snapshots[0][index], initial, abs_tol=1e-12):
                raise AssertionError(f"Incorrect initial condition at node {index}")
            if not math.isclose(snapshots[-1][index], final_value, rel_tol=1e-9, abs_tol=1e-11):
                raise AssertionError(f"Restarted time-series solution differs at node {index}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: verify_vtk_series_3d.py <nova_transient_heat_3d_executable>")
    verify(Path(sys.argv[1]).resolve())
    print("3D PVD/VTK frame sequence and monolithic final state verified.")
