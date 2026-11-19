#!/usr/bin/env python3
"""Independent end-to-end validation of streamed 3D PVD/VTK thermal output.

Uses only the Python standard library, with outputs in a temporary directory.
"""
from __future__ import annotations

import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def run(executable: Path, *args: str) -> str:
    result = subprocess.run(
        [str(executable), *args],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def read_vtk_scalar(path: Path) -> tuple[list[tuple[float, float, float]], list[float]]:
    lines = path.read_text(encoding="ascii").splitlines()
    assert "DATASET UNSTRUCTURED_GRID" in lines, path
    header = next(line for line in lines if line.startswith("POINTS "))
    count = int(header.split()[1])
    start = lines.index(header) + 1
    points = [tuple(map(float, line.split())) for line in lines[start:start + count]]
    assert all(len(p) == 3 and all(map(math.isfinite, p)) for p in points)
    assert f"POINT_DATA {count}" in lines
    scalar_header = lines.index("SCALARS temperature double 1")
    assert lines[scalar_header + 1] == "LOOKUP_TABLE default"
    field = [float(line) for line in lines[scalar_header + 2:]]
    assert len(field) == count
    assert all(map(math.isfinite, field))
    return points, field


def temperature_range(stdout: str) -> tuple[float, float]:
    match = re.search(
        r"Temperature range \[K\]: ([0-9.eE+-]+) \.\.\. ([0-9.eE+-]+)",
        stdout,
    )
    if not match:
        raise AssertionError("Missing documented final temperature range")
    return float(match.group(1)), float(match.group(2))


def verify(executable: Path) -> None:
    baseline = run(executable)
    with tempfile.TemporaryDirectory(prefix="nova-thermal-pvd-") as directory:
        target = Path(directory) / "thermal_case.pvd"
        result = run(executable, "--vtk-series", str(target))
        assert "ParaView collection:" in result
        assert temperature_range(result) == temperature_range(baseline)

        root = ET.parse(target).getroot()
        assert root.tag == "VTKFile" and root.get("type") == "Collection"
        datasets = root.findall("./Collection/DataSet")
        assert len(datasets) == 6  # initial plus five solved steps

        original_nodes = None
        for index, dataset in enumerate(datasets):
            filename = f"thermal_case_{index:06d}.vtk"
            assert dataset.get("file") == filename
            assert dataset.get("group") == ""
            actual_time = float(dataset.attrib["timestep"])
            assert math.isclose(actual_time, 20.0 * index, abs_tol=1e-12)
            points, values = read_vtk_scalar(target.parent / filename)
            if original_nodes is None:
                original_nodes = points
            else:
                assert points == original_nodes
            assert len(values) == len(points)
            for node, value in zip(points, values):
                if node[0] == 0.0:  # fixed-temperature wall
                    assert math.isclose(value, 300.0, abs_tol=1e-10)
            if index == 0:
                assert all(math.isclose(v, 300.0, abs_tol=1e-12) for v in values)
            if index == 5:
                lo, hi = temperature_range(baseline)
                assert math.isclose(min(values), lo, abs_tol=2e-8)
                assert math.isclose(max(values), hi, abs_tol=2e-8)

        invalid = subprocess.run(
            [str(executable), "--vtk-series", str(target.with_suffix(".txt"))],
            capture_output=True, text=True,
        )
        assert invalid.returncode != 0, "Invalid series extension was accepted"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: verify_transient_heat_3d_pvd.py <executable>")
    verify(Path(sys.argv[1]).resolve())
    print("Verified all 3D transient PVD/VTK frames against monolithic solver output.")
