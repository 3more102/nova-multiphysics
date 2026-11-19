#!/usr/bin/env python3
"""Independent end-to-end verification of tagged Gmsh transient PVD output."""

import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


def invoke(executable, mesh, *args, good=True):
    result = subprocess.run(
        [str(executable), str(mesh), *map(str, args)],
        capture_output=True,
        text=True,
        check=False,
    )
    if good and result.returncode:
        raise AssertionError("CLI failed:\n" + result.stderr + "\n" + result.stdout)
    if not good and result.returncode == 0:
        raise AssertionError("CLI accepted invalid arguments: " + repr(args))
    return result.stdout


def report_number(output, key):
    match = re.search(r"\b" + re.escape(key) + r"\s+([-+0-9.eE]+)", output)
    if match is None:
        raise AssertionError("Missing " + key + " in CLI report:\n" + output)
    value = float(match.group(1))
    if not math.isfinite(value):
        raise AssertionError("Non-finite " + key)
    return value


def read_vtk(path):
    lines = path.read_text(encoding="ascii").splitlines()
    if "DATASET UNSTRUCTURED_GRID" not in lines:
        raise AssertionError("Not an unstructured-grid VTK file: " + str(path))

    points_header = next(i for i, line in enumerate(lines) if line.startswith("POINTS "))
    point_count = int(lines[points_header].split()[1])
    points = [tuple(map(float, line.split()))
              for line in lines[points_header + 1:points_header + 1 + point_count]]
    if len(points) != point_count or any(
        len(point) != 3 or not all(map(math.isfinite, point)) for point in points
    ):
        raise AssertionError("Malformed VTK point coordinates")

    cells_header = next(i for i, line in enumerate(lines) if line.startswith("CELLS "))
    cell_count = int(lines[cells_header].split()[1])
    cells = []
    for line in lines[cells_header + 1:cells_header + 1 + cell_count]:
        values = list(map(int, line.split()))
        if len(values) != 5 or values[0] != 4:
            raise AssertionError("Expected four-node tetrahedral connectivity")
        cells.append(tuple(values[1:]))
    if len(cells) != cell_count:
        raise AssertionError("Malformed VTK cell connectivity")

    point_data = lines.index("POINT_DATA " + str(point_count))
    temperature_header = lines.index("SCALARS temperature double 1", point_data)
    temperature_lookup = lines.index("LOOKUP_TABLE default", temperature_header)
    temperatures = [float(value) for value in
                    lines[temperature_lookup + 1:temperature_lookup + 1 + point_count]]

    cell_data = lines.index("CELL_DATA " + str(cell_count))
    material_header = lines.index("SCALARS material_id int 1", cell_data)
    material_lookup = lines.index("LOOKUP_TABLE default", material_header)
    materials = [int(value) for value in
                 lines[material_lookup + 1:material_lookup + 1 + cell_count]]

    generation_header = lines.index("SCALARS heat_generation double 1", cell_data)
    generation_lookup = lines.index("LOOKUP_TABLE default", generation_header)
    generation = [float(value) for value in
                  lines[generation_lookup + 1:generation_lookup + 1 + cell_count]]

    flux_header = lines.index("VECTORS heat_flux double", cell_data)
    flux = [tuple(map(float, line.split()))
            for line in lines[flux_header + 1:flux_header + 1 + cell_count]]

    if (len(temperatures) != point_count or len(materials) != cell_count or
            len(generation) != cell_count or len(flux) != cell_count):
        raise AssertionError("VTK scalar/vector field length mismatch")
    numeric = temperatures + generation + [component for vector in flux for component in vector]
    if not all(map(math.isfinite, numeric)):
        raise AssertionError("VTK contains non-finite physical fields")
    if any(len(vector) != 3 for vector in flux):
        raise AssertionError("Malformed heat-flux vector")
    return points, cells, temperatures, materials, generation, flux


def read_series(path, expected_times):
    root = ET.parse(path).getroot()
    if root.tag != "VTKFile" or root.get("type") != "Collection":
        raise AssertionError("PVD root is not a VTK Collection")
    datasets = root.findall("./Collection/DataSet")
    if len(datasets) != len(expected_times):
        raise AssertionError("PVD frame count mismatch")
    frames = []
    reference_points = None
    reference_cells = None
    prefix = path.stem
    for index, (dataset, expected_time) in enumerate(zip(datasets, expected_times)):
        filename = f"{prefix}_{index:06d}.vtk"
        if dataset.get("file") != filename or Path(filename).name != filename:
            raise AssertionError("PVD contains an unexpected or unsafe frame name")
        if dataset.get("group") != "" or dataset.get("part") != "0":
            raise AssertionError("Unexpected PVD dataset metadata")
        actual_time = float(dataset.attrib["timestep"])
        if not math.isclose(actual_time, expected_time, rel_tol=0.0, abs_tol=1e-12):
            raise AssertionError("PVD physical timestep mismatch")
        frame = read_vtk(path.parent / filename)
        if reference_points is None:
            reference_points, reference_cells = frame[0], frame[1]
        elif frame[0] != reference_points or frame[1] != reference_cells:
            raise AssertionError("Mesh topology changed between time frames")
        frames.append(frame)
    return frames


def close_vectors(left, right, tolerance=1e-11):
    if len(left) != len(right):
        return False
    return all(abs(a - b) <= tolerance for a, b in zip(left, right))


def main():
    executable, wall_mesh, named_mesh = map(Path, sys.argv[1:4])
    with tempfile.TemporaryDirectory(prefix="nova-gmsh-transient-pvd-") as directory:
        root = Path(directory)

        # Exact uniform source ramp on two physical volume groups.
        uniform_pvd = root / "uniform.pvd"
        uniform_args = (
            "--material", "11", "1", "--material", "22", "3",
            "--capacity", "11", "2", "--capacity", "22", "2",
            "--initial", "300", "--start-time", "2", "--source", "4",
            "--dt", "0.1", "--steps", "5", "--scheme", "cn",
        )
        output = invoke(executable, wall_mesh, *uniform_args,
                        "--vtk-series", uniform_pvd)
        if abs(report_number(output, "final_time") - 2.5) > 1e-12:
            raise AssertionError("Wrong final physical time")
        if report_number(output, "max_abs_energy_balance_residual") > 1e-7:
            raise AssertionError("Uniform heating energy budget did not close")
        uniform_frames = read_series(uniform_pvd, [2.0 + 0.1 * i for i in range(6)])
        for index, frame in enumerate(uniform_frames):
            expected = 300.0 + 0.2 * index
            if any(abs(value - expected) > 1e-7 for value in frame[2]):
                raise AssertionError("Analytical uniform temperature disagrees with VTK frame")
            if set(frame[3]) != {11, 22}:
                raise AssertionError("VTK lost one or more Gmsh volume tags")
            if any(abs(value - 4.0) > 1e-12 for value in frame[4]):
                raise AssertionError("VTK heat-generation field is incorrect")
            if any(math.sqrt(sum(component * component for component in vector)) > 1e-9
                   for vector in frame[5]):
                raise AssertionError("Uniform temperature should have zero conductive flux")

        # The final PVD frame must equal the independently exported final VTK state.
        final_vtk = root / "uniform_final.vtk"
        invoke(executable, wall_mesh, *uniform_args, "--vtk", final_vtk)
        final_frame = read_vtk(final_vtk)
        if not close_vectors(uniform_frames[-1][2], final_frame[2]):
            raise AssertionError("PVD final temperature differs from final-only VTK export")
        if uniform_frames[-1][3:] != final_frame[3:]:
            raise AssertionError("PVD final cell fields differ from final-only VTK export")

        # Mixed scalar/tensor regions, regional heating, Dirichlet projection and convection.
        mixed_pvd = root / "mixed.pvd"
        mixed_args = (
            "--material", "11", "1",
            "--material-tensor", "22", "3", "2", "1", "0", "0", "0",
            "--capacity", "11", "2", "--capacity", "22", "4",
            "--initial", "300", "--start-time", "1",
            "--temperature", "101", "295", "--convection", "102", "10", "290",
            "--heat-source", "22", "0.5", "--dt", "0.05", "--steps", "3",
            "--scheme", "be", "--vtk-series", mixed_pvd,
        )
        mixed_output = invoke(executable, wall_mesh, *mixed_args)
        if report_number(mixed_output, "max_abs_energy_balance_residual") > 1e-7:
            raise AssertionError("Mixed-material energy budget did not close")
        mixed_frames = read_series(mixed_pvd, [1.0, 1.05, 1.10, 1.15])
        for frame in mixed_frames:
            if min(frame[2]) < 289.9 or max(frame[2]) > 301.0:
                raise AssertionError("Mixed-material temperatures outside expected bounds")
            if set(frame[3]) != {11, 22} or set(frame[4]) != {0.0, 0.5}:
                raise AssertionError("Mixed-material VTK cell fields are incorrect")
        for point, value in zip(mixed_frames[0][0], mixed_frames[0][2]):
            if abs(point[0]) < 1e-12 and abs(value - 295.0) > 1e-10:
                raise AssertionError("Initial frame omitted Dirichlet projection")

        # Resolve Gmsh physical names as well as numeric tags.
        named_pvd = root / "named.pvd"
        named_args = (
            "--material", "solid", "1", "--capacity", "solid", "2",
            "--initial", "300", "--temperature", "fixed-temperature", "310",
            "--dt", "0.25", "--steps", "2", "--scheme", "be",
            "--vtk-series", named_pvd,
        )
        named_output = invoke(executable, named_mesh, *named_args)
        if report_number(named_output, "max_abs_energy_balance_residual") > 1e-7:
            raise AssertionError("Named-group energy budget did not close")
        named_frames = read_series(named_pvd, [0.0, 0.25, 0.5])
        if any(tag != 42 for tag in named_frames[0][3]):
            raise AssertionError("Named volume group did not reach VTK")
        for point, value in zip(named_frames[0][0], named_frames[0][2]):
            if abs(point[2]) < 1e-12 and abs(value - 310.0) > 1e-10:
                raise AssertionError("Named Dirichlet group was not projected at t=0")

        # Malformed requests must fail before publishing a collection.
        invoke(executable, wall_mesh, *uniform_args, "--vtk-series",
               root / "wrong_extension.txt", good=False)
        invoke(executable, wall_mesh, *uniform_args, "--vtk-series",
               root / "conflict.pvd", "--vtk", root / "conflict.vtk", good=False)
        invoke(executable, wall_mesh,
               "--material", "11", "1", "--capacity", "11", "2",
               "--initial", "300", "--dt", "0.1", "--steps", "1", good=False)
        bad_capacity_args = (
            "--material", "11", "1", "--material", "22", "3",
            "--capacity", "11", "-1", "--capacity", "22", "2",
            "--initial", "300", "--dt", "0.1", "--steps", "1",
        )
        invoke(executable, wall_mesh, *bad_capacity_args, good=False)
        invoke(executable, wall_mesh, *uniform_args, "--temperature", "999", "300", good=False)
        bad_tensor_args = (
            "--material", "11", "1",
            "--material-tensor", "22", "-1", "1", "1", "0", "0", "0",
            "--capacity", "11", "2", "--capacity", "22", "2",
            "--initial", "300", "--dt", "0.1", "--steps", "1",
        )
        failed_pvd = root / "invalid_tensor.pvd"
        invoke(executable, wall_mesh, *bad_tensor_args, "--vtk-series", failed_pvd, good=False)
        if failed_pvd.exists() or list(root.glob(".invalid_tensor.nova-tmp-*")):
            raise AssertionError("Failed solve left a partial PVD series behind")
        invoke(executable, wall_mesh, *uniform_args, "--vtk-series", uniform_pvd, good=False)

    print("PASS: Gmsh transient PVD frames, numerical fields, physical times, names and invalid inputs")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: verify_gmsh_transient_heat_3d_pvd.py <executable> <wall.msh> <named.msh>")
    main()
