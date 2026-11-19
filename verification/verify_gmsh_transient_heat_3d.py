#!/usr/bin/env python3
"""End-to-end verification of imported, tagged transient 3D thermal studies."""

import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def invoke(exe, mesh, *args, good=True):
    run = subprocess.run([exe, str(mesh), *map(str, args)],
                         capture_output=True, text=True, check=False)
    if good and run.returncode:
        raise AssertionError("CLI failed: " + run.stderr + "\n" + run.stdout)
    if not good and run.returncode == 0:
        raise AssertionError("CLI accepted invalid inputs: " + str(args))
    return run.stdout


def report_number(text, key):
    match = re.search(r"\b" + re.escape(key) +
                      r"\s+([-+0-9.eE]+)", text)
    if match is None:
        raise AssertionError("Missing " + key + " from report:\n" + text)
    value = float(match.group(1))
    if not math.isfinite(value):
        raise AssertionError("Non-finite " + key)
    return value


def temperatures(path):
    vtk = path.read_text(encoding="utf-8")
    match = re.search(r"POINT_DATA\s+(\d+)\s+SCALARS\s+temperature\s+"
                      r"\w+(?:\s+\d+)?\s+LOOKUP_TABLE\s+default\s+",
                      vtk, re.IGNORECASE)
    if match is None:
        raise AssertionError("VTK missing nodal temperature field")
    count = int(match.group(1))
    fields = vtk[match.end():].split()
    values = [float(item) for item in fields[:count]]
    if len(values) != count or not all(map(math.isfinite, values)):
        raise AssertionError("Malformed VTK temperatures")
    if "CELL_DATA" not in vtk or "heat_generation" not in vtk:
        raise AssertionError("VTK missing material/cell heating data")
    return values


def main():
    exe, mesh = sys.argv[1:3]
    with tempfile.TemporaryDirectory(prefix="nova-gmsh-transient-") as temp:
        vtk = Path(temp) / "uniform_heat.vtk"
        args = ("--material", 11, 1, "--material", 22, 3,
                "--capacity", 11, 2, "--capacity", 22, 2,
                "--initial", 300, "--source", 4, "--dt", 0.1,
                "--steps", 5, "--scheme", "cn")
        uniform = invoke(exe, mesh, *args, "--vtk", vtk)
        assert abs(report_number(uniform, "final_time") - 0.5) < 1e-12
        assert report_number(uniform, "max_abs_energy_balance_residual") < 1e-7
        values = temperatures(vtk)
        if max(abs(value - 301.0) for value in values) > 1e-7:
            raise AssertionError("Analytical uniform source ramp disagrees with VTK")

        # Test two tagged exterior surfaces and a mixed material model.
        mixed = invoke(
            exe, mesh, "--material", 11, 1,
            "--material-tensor", 22, 3, 2, 1, 0, 0, 0,
            "--capacity", 11, 2, "--capacity", 22, 4,
            "--initial", 300, "--temperature", 101, 295,
            "--convection", 102, 10, 290,
            "--dt", 0.05, "--steps", 3, "--scheme", "be"
        )
        assert report_number(mixed, "min_temperature") >= 289.9
        assert report_number(mixed, "max_temperature") <= 300.1
        assert report_number(mixed, "max_abs_energy_balance_residual") < 1e-7

        # Missing volume coverage and invalid physical parameters must fail.
        invoke(exe, mesh, "--material", 11, 1, "--capacity", 11, 2,
               "--capacity", 22, 2, "--initial", 300, "--dt", 0.1,
               "--steps", 1, good=False)
        invoke(exe, mesh, *args, "--capacity", 11, -1, good=False)
        invoke(exe, mesh, *args, "--temperature", 999, 300, good=False)
        invoke(exe, mesh, *args, "--vtk", mesh, good=False)

    print("PASS: tagged Gmsh heat ramp, tensor convection, balances, VTK, invalid inputs")


if __name__ == "__main__":
    main()
