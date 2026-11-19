#!/usr/bin/env python3
"""Independent end-to-end verification of single-assembly 3D thermal PVD output."""

import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def scalar_values(path: Path) -> list[float]:
    lines = path.read_text(encoding="ascii").splitlines()
    point_header = next(line for line in lines if line.startswith("POINTS "))
    count = int(point_header.split()[1])
    marker = lines.index("SCALARS temperature double 1")
    if lines[marker + 1] != "LOOKUP_TABLE default":
        raise AssertionError("Missing VTK scalar lookup table")
    values = [float(value) for value in lines[marker + 2:marker + 2 + count]]
    if len(values) != count or not all(math.isfinite(v) for v in values):
        raise AssertionError("Invalid nodal temperatures")
    return values


def main(executable: str) -> None:
    with tempfile.TemporaryDirectory(prefix="nova_thermal_pvd_") as directory:
        folder = Path(directory)
        process = subprocess.run(
            [executable, "--vtk-series", "thermal_verify"],
            cwd=folder,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if "PVD collection written: thermal_verify.pvd (6 frames)" not in process.stdout:
            raise AssertionError("Expected exactly six snapshots including initial state")
        root = ET.parse(folder / "thermal_verify.pvd").getroot()
        entries = root.find("Collection")
        if root.tag != "VTKFile" or entries is None:
            raise AssertionError("Not a VTK time collection")
        frames = entries.findall("DataSet")
        if len(frames) != 6:
            raise AssertionError(f"Expected 6 frames, got {len(frames)}")
        all_values = []
        for step, entry in enumerate(frames):
            time = float(entry.attrib["timestep"])
            if not math.isclose(time, 20.0 * step, abs_tol=1e-12):
                raise AssertionError(f"Incorrect timestep {time} at frame {step}")
            name = entry.attrib["file"]
            if name != f"thermal_verify_{step:06d}.vtk":
                raise AssertionError("Unsafe or incorrectly ordered VTK frame reference")
            all_values.append(scalar_values(folder / name))
        if len({len(frame) for frame in all_values}) != 1:
            raise AssertionError("Changing mesh vertex count")
        if any(abs(temp - 300.0) > 1e-10 for temp in all_values[0]):
            raise AssertionError("Initial exported state is not 300 K")
        if max(all_values[-1]) <= 300.0:
            raise AssertionError("Heated solid did not warm")
        match = re.search(
            r"Temperature range \[K\]: ([-+0-9.eE]+) \.\.\. ([-+0-9.eE]+)",
            process.stdout,
        )
        if match is None:
            raise AssertionError("Missing final solver temperature summary")
        low, high = (float(item) for item in match.groups())
        if not math.isclose(min(all_values[-1]), low, abs_tol=1e-7) or not math.isclose(
            max(all_values[-1]), high, abs_tol=1e-7
        ):
            raise AssertionError("Final exported snapshot disagrees with solver result")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: verify_transient_series_3d.py /path/to/nova_transient_heat_3d")
    main(sys.argv[1])
