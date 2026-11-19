#!/usr/bin/env python3
"""Verify current transient tetrahedral example's ParaView PVD and VTK frames."""
from __future__ import annotations

import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def read_vtk_temperature(path: Path) -> list[float]:
    lines = path.read_text(encoding="ascii").splitlines()
    if "DATASET UNSTRUCTURED_GRID" not in lines:
        raise AssertionError("Not a VTK unstructured grid")
    n = int(next(line.split()[1] for line in lines if line.startswith("POINTS ")))
    if f"POINT_DATA {n}" not in lines:
        raise AssertionError("Bad VTK point data size")
    loc = lines.index("LOOKUP_TABLE default")
    values = [float(value) for value in lines[loc + 1:]]
    if len(values) != n or not all(math.isfinite(x) for x in values):
        raise AssertionError("Bad VTK temperature field")
    return values


def verify(exe: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="nova-thermal3d-pvd-") as temp:
        folder = Path(temp)
        p = subprocess.run([str(exe), "--vtk-series", "cooling"], cwd=folder,
                           text=True, capture_output=True, check=True)
        collection = folder / "cooling.pvd"
        root = ET.parse(collection).getroot()
        if root.tag != "VTKFile" or root.attrib.get("type") != "Collection":
            raise AssertionError("PVD must be a VTK Collection")
        entries = root.findall("./Collection/DataSet")
        if len(entries) != 6:
            raise AssertionError("Expected initial and five thermal frames")
        frames = []
        for i, entry in enumerate(entries):
            expected = f"cooling_{i:06d}.vtk"
            if entry.attrib.get("file") != expected:
                raise AssertionError(f"Incorrect frame filename {i}")
            if not math.isclose(float(entry.attrib["timestep"]), i * 20.0,
                                abs_tol=1e-12):
                raise AssertionError(f"Incorrect time for frame {i}")
            frames.append(read_vtk_temperature(folder / expected))
        if not all(abs(x - 300.0) <= 1e-12 for x in frames[0]):
            raise AssertionError("Initial temperature should be 300 K")
        if any(len(x) != len(frames[0]) for x in frames):
            raise AssertionError("Inconsistent mesh vertex count across frames")
        m = re.search(r"Temperature range \\[K\\]:\\s*([0-9.eE+-]+)\\s*\\.\\.\\.\\s*([0-9.eE+-]+)", p.stdout)
        if not m:
            raise AssertionError("Missing solver temperature range")
        lo, hi = map(float, m.groups())
        if not (math.isclose(min(frames[-1]), lo, abs_tol=1e-8) and
                math.isclose(max(frames[-1]), hi, abs_tol=1e-8)):
            raise AssertionError("Last ParaView frame differs from solver result")
        invalid = subprocess.run([str(exe), "--vtk-series", "../unsafe"], cwd=folder,
                                 text=True, capture_output=True)
        if invalid.returncode == 0 or (folder.parent / "unsafe.pvd").exists():
            raise AssertionError("Unsafe output prefix was accepted")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: verify_transient_pvd_3d.py executable")
    verify(Path(sys.argv[1]).resolve())
    print("Current 3D transient VTK/PVD time-series verified.")
