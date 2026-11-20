#!/usr/bin/env python3
"""End-to-end verification of imported, tagged transient 3D thermal studies."""

import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def invoke(exe, mesh, *args, good=True):
    run = subprocess.run([str(exe), str(mesh), *map(str, args)],
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
    if "CELL_DATA" not in vtk or "material_id" not in vtk:
        raise AssertionError("VTK missing cell material data")
    if "heat_generation" not in vtk or "heat_flux" not in vtk:
        raise AssertionError("VTK missing cell heating or heat-flux data")
    return values


def verify_series(exe, mesh, args, pvd):
    report = invoke(exe, mesh, *args, "--vtk-series", pvd)
    assert abs(report_number(report, "final_time") - 0.5) < 1e-12
    assert report_number(report, "max_abs_energy_balance_residual") < 1e-7

    root = ET.parse(pvd).getroot()
    assert root.tag == "VTKFile" and root.get("type") == "Collection"
    datasets = root.findall("./Collection/DataSet")
    assert len(datasets) == 6  # projected initial state plus five solved steps

    original_names = []
    original_points = None
    for index, dataset in enumerate(datasets):
        expected_name = f"uniform_heat_{index:06d}.vtk"
        assert dataset.get("file") == expected_name
        assert dataset.get("group") == ""
        time = float(dataset.attrib["timestep"])
        expected_time = 0.1 * index
        assert math.isclose(time, expected_time, rel_tol=0.0, abs_tol=1e-14)
        frame = pvd.parent / expected_name
        assert frame.is_file()
        text = frame.read_text(encoding="utf-8")
        point_match = re.search(r"POINTS\s+(\d+)\s+\w+\s+", text)
        if point_match is None:
            raise AssertionError("VTK frame has no point geometry")
        point_count = int(point_match.group(1))
        point_tokens = text[point_match.end():].split()
        points = tuple(tuple(map(float, point_tokens[i:i + 3]))
                       for i in range(0, 3 * point_count, 3))
        if len(points) != point_count or not all(
            len(point) == 3 and all(map(math.isfinite, point)) for point in points
        ):
            raise AssertionError("Malformed VTK point geometry")
        if original_points is None:
            original_points = points
        elif points != original_points:
            raise AssertionError("Mesh geometry changed between transient frames")

        values = temperatures(frame)
        expected_temperature = 300.0 + 2.0 * expected_time
        if max(abs(value - expected_temperature) for value in values) > 1e-7:
            raise AssertionError(
                f"Frame {index} disagrees with analytical T(t)=300+2t"
            )
        original_names.append(expected_name)

    # A second run must not truncate an existing collection or any of its frames.
    invoke(exe, mesh, *args, "--vtk-series", pvd, good=False)
    assert ET.parse(pvd).getroot().tag == "VTKFile"
    assert all((pvd.parent / name).is_file() for name in original_names)


def verify_mesh_quality(exe, mesh, thermal_args, directory):
    # Inspection must work without thermal parameters or a solver run.
    report = invoke(exe, mesh, "--check-mesh")
    ratio = report_number(report, "Minimum mean ratio:")
    volume = report_number(report, "Total tetrahedral volume [m^3]:")
    angle = report_number(report, "Minimum dihedral angle [degrees]:")
    if not (0.0 < ratio < 1.0 and volume > 0.0 and angle > 0.0):
        raise AssertionError("Invalid tetrahedral quality report")
    if "Worst tetrahedron index (zero-based):" not in report:
        raise AssertionError("Missing worst-cell diagnostic")

    accepted = format(0.5 * ratio, ".17g")
    rejected = format(0.5 * (1.0 + ratio), ".17g")
    assert "Minimum mean ratio:" in invoke(
        exe, mesh, "--min-mean-ratio", accepted, "--check-mesh")

    # Report/gate option order must not affect an otherwise valid heat solve.
    for quality_options in (
        ("--quality-report", "--min-mean-ratio", accepted),
        ("--min-mean-ratio", accepted, "--quality-report"),
    ):
        solved = invoke(exe, mesh, *thermal_args, *quality_options)
        assert "Minimum mean ratio:" in solved
        assert "final_time" in solved

    # Reject the mesh before writing any ParaView frame or collection.
    rejected_pvd = directory / "blocked_quality.pvd"
    command = [str(exe), str(mesh), *map(str, thermal_args),
               "--min-mean-ratio", rejected, "--vtk-series", str(rejected_pvd)]
    failed = subprocess.run(command, capture_output=True, text=True, check=False)
    if failed.returncode == 0 or "quality gate failed" not in failed.stderr:
        raise AssertionError("Quality gate failed to reject the transient solve")
    if rejected_pvd.exists() or list(directory.glob("blocked_quality_*.vtk")):
        raise AssertionError("Quality-gate failure created partial result files")

    for invalid in ("0", "-0.1", "1.001", "nan", "not_a_number"):
        invoke(exe, mesh, "--check-mesh", "--min-mean-ratio", invalid,
               good=False)
    invoke(exe, mesh, "--check-mesh", "--min-mean-ratio", accepted,
           "--min-mean-ratio", accepted, good=False)
    invoke(exe, mesh, "--check-mesh", "--quality-report",
           "--quality-report", good=False)
    invoke(exe, mesh, "--check-mesh", "--material", 11, 1, good=False)
    invoke(exe, mesh, "--check-mesh", "--scheme", "cn", good=False)


def main():
    exe, mesh = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory(prefix="nova-gmsh-transient-") as temp:
        directory = Path(temp)
        args = ("--material", 11, 1, "--material", 22, 3,
                "--capacity", 11, 2, "--capacity", 22, 2,
                "--initial", 300, "--source", 4, "--dt", 0.1,
                "--steps", 5, "--scheme", "cn")
        verify_mesh_quality(exe, mesh, args, directory)
        pvd = directory / "uniform_heat.pvd"
        verify_series(exe, mesh, args, pvd)

        final_vtk = directory / "uniform_final.vtk"
        final_report = invoke(exe, mesh, *args, "--vtk", final_vtk)
        assert abs(report_number(final_report, "final_time") - 0.5) < 1e-12
        final_values = temperatures(final_vtk)
        assert max(abs(value - 301.0) for value in final_values) < 1e-7

        # No source and no boundary exchange preserve a uniform insulated field.
        equilibrium = invoke(
            exe, mesh, "--material", 11, 1, "--material", 22, 3,
            "--capacity", 11, 2, "--capacity", 22, 2,
            "--initial", 300, "--dt", 0.1, "--steps", 3
        )
        assert math.isclose(report_number(equilibrium, "min_temperature"),
                            300.0, abs_tol=1e-9)
        assert math.isclose(report_number(equilibrium, "max_temperature"),
                            300.0, abs_tol=1e-9)
        assert report_number(equilibrium, "max_abs_energy_balance_residual") < 1e-8

        # A cell-only source is additive and must increase total thermal energy.
        cell_source = invoke(
            exe, mesh, "--material", 11, 1, "--material", 22, 3,
            "--capacity", 11, 2, "--capacity", 22, 2,
            "--initial", 300, "--heat-source", 22, 4,
            "--dt", 0.1, "--steps", 5
        )
        assert report_number(cell_source, "max_temperature") > 300.0
        assert report_number(cell_source, "max_abs_energy_balance_residual") < 1e-7

        # Positive prescribed outward flux cools the tagged exterior face.
        neumann = invoke(
            exe, mesh, "--material", 11, 1, "--material", 22, 3,
            "--capacity", 11, 2, "--capacity", 22, 2,
            "--initial", 300, "--flux", 102, 2,
            "--dt", 0.05, "--steps", 3, "--scheme", "be"
        )
        assert report_number(neumann, "min_temperature") < 300.0
        report_number(neumann, "max_temperature")  # finite, even if the P1 mesh overshoots locally
        assert report_number(neumann, "max_abs_energy_balance_residual") < 1e-7

        # Mixed scalar/tensor materials, fixed temperature and convection.
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

        # Missing material coverage, invalid data, tags, and output modes fail.
        invoke(exe, mesh, "--material", 11, 1,
               "--capacity", 11, 2, "--capacity", 22, 2,
               "--initial", 300, "--dt", 0.1, "--steps", 1, good=False)
        bad_capacity = (
            "--material", 11, 1, "--material", 22, 3,
            "--capacity", 11, -1, "--capacity", 22, 2,
            "--initial", 300, "--dt", 0.1, "--steps", 1,
        )
        invoke(exe, mesh, *bad_capacity, good=False)
        invoke(exe, mesh, *args, "--temperature", 999, 300, good=False)
        invoke(exe, mesh, *args, "--temperature", 101, 300,
               "--flux", 101, 1, good=False)
        bad_tensor = (
            "--material", 11, 1,
            "--material-tensor", 22, 1, 1, 1, 2, 0, 0,
            "--capacity", 11, 2, "--capacity", 22, 2,
            "--initial", 300, "--dt", 0.1, "--steps", 1,
        )
        invoke(exe, mesh, *bad_tensor, good=False)
        invoke(exe, mesh, *args, "--vtk", mesh, good=False)
        invoke(exe, mesh, *args, "--vtk-series",
               directory / "invalid.txt", good=False)
        invoke(exe, mesh, *args, "--vtk", directory / "both.vtk",
               "--vtk-series", directory / "both.pvd", good=False)

    print("PASS: Gmsh transient BE/CN, quality gates, PVD frames, cell fields, source, flux, energy and invalid-input checks")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: verify_gmsh_transient_heat_3d.py <executable> <mesh.msh>")
    main()
