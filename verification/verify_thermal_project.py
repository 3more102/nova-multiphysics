"""Independent, standard-library checks of the installed nova-run file protocol."""
import argparse
import math
from pathlib import Path
import subprocess
import tempfile


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def near(actual, reference, tolerance, description):
    require(math.isfinite(actual) and abs(actual - reference) <= tolerance,
            f"{description}: {actual:.12g}, reference {reference:.12g}, tolerance {tolerance:g}")


def invoke(runner, *args, success=True):
    cmd = [str(runner), *map(str, args)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
    if success:
        require(result.returncode == 0, f"{cmd}: {result.stderr}")
    else:
        require(result.returncode > 0 and result.stderr and not result.stdout,
                f"Unexpected failure response: {result}")
    return result.stdout


def read_vtk(path):
    tokens = path.read_text(encoding="ascii").split()

    def numbers(key, offset, width, count, cast=float):
        i = tokens.index(key) + offset
        return [tuple(map(cast, tokens[j:j + width])) for j in range(i, i + count * width, width)]

    count = int(tokens[tokens.index("POINTS") + 1])
    cells = int(tokens[tokens.index("CELLS") + 1])
    points = numbers("POINTS", 3, 3, count)
    connectivity = numbers("CELLS", 3, 5, cells, int)
    require(all(c[0] == 4 and len(set(c[1:])) == 4 and all(0 <= n < count for n in c[1:])
                for c in connectivity), "Invalid tetrahedron in exported mesh")
    require(int(tokens[tokens.index("POINT_DATA") + 1]) == count, "VTK point count")
    require(int(tokens[tokens.index("CELL_DATA") + 1]) == cells, "VTK cell count")
    temperatures = [x[0] for x in numbers("temperature", 5, 1, count)]
    fluxes = numbers("heat_flux", 2, 3, cells)
    return points, [c[1:] for c in connectivity], temperatures, fluxes


def solve(runner, project, vtk):
    text = invoke(runner, "solve", project, "--vtk", vtk)
    lines = text.splitlines()
    require(lines[0] == "NOVA_RESULT 1", "Unknown result protocol")
    fields = dict(line.split() for line in lines[1:])
    data = read_vtk(vtk)
    require(len(data[0]) == int(fields["nodes"]) and len(data[1]) == int(fields["tetrahedra"]),
            "Process summary and VTK disagree")
    return {k: float(v) for k, v in fields.items() if k not in {"project", "study"}}, data


def tetrahedron_volume(points, cell):
    """Geometric volume of a linear tetrahedron, independent of node winding."""
    p0, p1, p2, p3 = (points[i] for i in cell)
    a = tuple(p1[j] - p0[j] for j in range(3))
    b = tuple(p2[j] - p0[j] for j in range(3))
    c = tuple(p3[j] - p0[j] for j in range(3))
    determinant = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    volume = abs(determinant) / 6.0
    require(math.isfinite(volume) and volume > 0.0, "Invalid tetrahedron volume in exported VTK")
    return volume


def wall(runner, fixture, work):
    # Thermal resistances per unit area in series. No NOVA or FEM computation.
    k, length, area, h, hot, ambient = 2.0, 0.5, 0.25, 12.0, 400.0, 300.0
    q = (hot - ambient) / (length / k + 1 / h)
    rate = area * q
    metrics, (points, _, temperature, fluxes) = solve(runner, fixture, work / "wall.vtk")
    error = max(abs(t - (hot - q * x / k)) for (x, _, _), t in zip(points, temperature))
    flux_error = max(math.dist(f, (q, 0, 0)) for f in fluxes)
    near(error, 0, 2e-8, "Maximum nodal field error [K]")
    near(flux_error, 0, 1e-6, "Maximum Fourier flux vector error [W/m^2]")
    near(metrics["convection_out_W"], rate, 2e-8, "Convective power [W]")
    near(metrics["fixed_out_W"], -rate, 2e-8, "Hot face reaction [W]")
    near(metrics["balance_defect_W"], 0, 2e-8, "Global balance [W]")
    near(metrics["min_K"], ambient + q / h, 2e-8, "Cold face [K]")
    print(f"wall: nodes={len(points)}, cold_K={metrics['min_K']:.14g}, convection_W={metrics['convection_out_W']:.14g}, "
          f"max_error_K={error:.6g}, max_flux_error_W_m2={flux_error:.6g}, balance_W={metrics['balance_defect_W']:.6g}")


def generated(runner, fixture, work):
    # -kT''=Q, T(0)=400 and -kT'(L)=g. The box has insulated lateral faces.
    length, k, source, g = 0.5, 2.0, 240.0, 40.0
    base = fixture.read_text(encoding="ascii").replace("source steady 0", "source steady 240").replace(
        "convection 12 300", "outward_flux 40")
    previous = None
    for n in (2, 4, 8, 16):
        path = work / "generation.nova"
        path.write_text(base.replace("8 3 2", f"{n} {n} {n}"), encoding="ascii")
        metrics, (points, cells, temperature, fluxes) = solve(runner, path, work / "generation.vtk")
        t_sq, q_sq = [], []
        for cell, flux in zip(cells, fluxes):
            x = sum(points[node][0] for node in cell) / 4
            computed = sum(temperature[node] for node in cell) / 4
            reference = 400 + ((source * length - g) * x - source * x * x / 2) / k
            t_sq.append((computed - reference) ** 2)
            q_sq.append(math.dist(flux, (g - source * (length - x), 0, 0)) ** 2)
        t_rms = math.sqrt(math.fsum(t_sq) / len(cells))
        q_rms = math.sqrt(math.fsum(q_sq) / len(cells))
        near(metrics["generated_W"], source * length**3, 2e-8, "Total generation [W]")
        near(metrics["neumann_out_W"], g * length**2, 2e-8, "Prescribed outward power [W]")
        near(metrics["fixed_out_W"], source * length**3 - g * length**2, 2e-8, "Reaction [W]")
        near(metrics["balance_defect_W"], 0, 2e-8, "Heat balance [W]")
        print(f"n={n} nodes={len(points)} tets={len(cells)} T_centroid_rms_K={t_rms:.8g} "
              f"q_centroid_rms_W_m2={q_rms:.8g} balance_W={metrics['balance_defect_W']:.6g}")
        if previous:
            t_order, q_order = math.log2(previous[0] / t_rms), math.log2(previous[1] / q_rms)
            require(t_order > 1.6 and q_order > 0.8, f"Refinement orders {t_order}, {q_order}")
            print(f"orders T={t_order:.5f} q={q_order:.5f}")
        previous = t_rms, q_rms


def cube_mean_reference(length, k, source, cutoff=127):
    """Positive volume-mean sine series; upper bound includes every omitted mode."""
    odd = range(1, cutoff + 1, 2)
    terms = (1.0 / (a*a * b*b * c*c * (a*a + b*b + c*c)) for a in odd for b in odd for c in odd)
    scale = source * length**2 / k
    lower = 512 * scale * math.fsum(terms) / math.pi**8
    tail = 8 * scale / (math.pi**4 * cutoff**3)
    return lower, lower + tail


def cube_fd_reference(n, length, k, source):
    """Independent seven-point FD solve using SciPy's orthonormal sine transform."""
    import numpy as np
    from scipy.fft import dstn, idstn

    require(n % 2 == 0, "The center probe requires an even FD subdivision count")
    h = length / n
    mode = np.arange(1, n, dtype=float)
    axis = 4 * np.sin(np.pi * mode / (2 * n))**2 / h**2
    eigenvalues = axis[:, None, None] + axis[None, :, None] + axis[None, None, :]
    rhs = np.full((n - 1, n - 1, n - 1), source / k)
    rise = idstn(dstn(rhs, type=1, norm="ortho") / eigenvalues, type=1, norm="ortho")
    # All boundary rises are zero in the tensor-product trapezoidal integral.
    mean_rise = float(np.sum(rise, dtype=float) / n**3)
    center_rise = float(rise[n // 2 - 1, n // 2 - 1, n // 2 - 1])
    return mean_rise, center_rise, rise


def cube(runner, fixture, work, check_scipy=False):
    length, k, source, boundary = 0.5, 2.0, 240.0, 300.0
    lower, upper = cube_mean_reference(length, k, source)
    print(f"Continuum mean rise bracket [{lower:.12g}, {upper:.12g}] K")
    base = fixture.read_text(encoding="ascii")
    mesh_record = "structured_tet 16 16 16"
    require(base.count(mesh_record) == 1, "Expected one 16-subdivision cube mesh record")
    previous_errors = None
    for n in (4, 8, 16):
        path = work / "cube.nova"
        path.write_text(base.replace(mesh_record, f"structured_tet {n} {n} {n}"), encoding="ascii")
        metrics, (points, cells, temperature, _) = solve(runner, path, work / "cube.vtk")
        boundary_error = max(abs(t - boundary) for point, t in zip(points, temperature)
                             if any(x in (0, length) for x in point))
        near(boundary_error, 0, 2e-8, "All fixed boundary nodes [K]")

        # The structured mesh is invariant under inversion through the cube center.
        # Check every nodal value against its geometrically reflected counterpart.
        point_index = {tuple(round(x, 12) for x in point): i for i, point in enumerate(points)}
        require(len(point_index) == len(points), "Cube mesh contains duplicate rounded coordinates")
        symmetry_error = 0.0
        for i, point in enumerate(points):
            reflected = tuple(round(length - x, 12) for x in point)
            partner = point_index.get(reflected)
            require(partner is not None, f"Missing reflected cube node for {point}")
            symmetry_error = max(symmetry_error, abs(temperature[i] - temperature[partner]))
        near(symmetry_error, 0, 2e-8, "Inversion symmetry of nodal temperature [K]")

        center_nodes = [i for i, point in enumerate(points)
                        if all(abs(x - length / 2) <= 1e-12 for x in point)]
        require(len(center_nodes) == 1, "Expected exactly one cube-center node")
        center_rise = temperature[center_nodes[0]] - boundary
        near(temperature[center_nodes[0]], max(temperature), 2e-8,
             "Maximum nodal temperature must occur at the cube center [K]")
        near(metrics["max_K"], temperature[center_nodes[0]], 2e-8,
             "Reported peak temperature must match the center probe [K]")

        # A linear P1 field integrates exactly as the mean of its four vertices.
        # Weight each element by its geometric volume; do not assume equal cells.
        volumes, weighted_rises = [], []
        for cell in cells:
            volume = tetrahedron_volume(points, cell)
            volumes.append(volume)
            weighted_rises.append(
                volume * math.fsum(temperature[i] - boundary for i in cell) / 4
            )
        total_volume = math.fsum(volumes)
        near(total_volume, length**3, 2e-12, "Exported tetrahedral cube volume [m^3]")
        mean_rise = math.fsum(weighted_rises) / total_volume
        errors = (lower - mean_rise, upper - mean_rise)
        require(errors[0] > 0, f"Expected refinement from below, observed {mean_rise:g} K")
        near(metrics["generated_W"], source * length**3, 2e-8, "Cube heat generation [W]")
        near(metrics["fixed_out_W"], source * length**3, 2e-8, "Six-face outward reaction [W]")
        near(metrics["balance_defect_W"], 0, 2e-8, "Cube heat balance [W]")
        print(f"cube n={n} mean_rise_K={mean_rise:.12g} center_rise_K={center_rise:.12g} "
              f"volume_m3={total_volume:.12g} symmetry_error_K={symmetry_error:.3g} "
              f"error_upper_K={errors[1]:.8g} "
              f"fixed_out_W={metrics['fixed_out_W']:.12g} balance_W={metrics['balance_defect_W']:.6g}")
        if previous_errors:
            # Worst-case order over both continuum error intervals.
            order_lower = math.log2(previous_errors[0] / errors[1])
            require(order_lower > 1.65, f"Cube mean convergence order lower bound {order_lower}")
            print(f"mean_order_lower={order_lower:.7g}")
        previous_errors = errors
    require(errors[1] < 0.025, "Finest cube mean exceeds its 0.025 K error allowance")
    if check_scipy:
        previous_mean = None
        previous_center = None
        previous_center_difference = None
        for n in (16, 32, 64):
            fd_mean, fd_center, fd_rise = cube_fd_reference(n, length, k, source)
            require(fd_mean < lower, "FD cube reference exceeds continuum bracket")
            if previous_mean is not None:
                require(math.log2((lower - previous_mean) / (upper - fd_mean)) > 1.8,
                        "FD cube mean refinement")
            if previous_center is not None:
                center_difference = abs(previous_center - fd_center)
                if previous_center_difference is not None:
                    center_order = math.log2(previous_center_difference / center_difference)
                    require(center_order > 1.8, f"FD cube center refinement order {center_order}")
                    print(f"center_refinement_order={center_order:.7g}")
                previous_center_difference = center_difference
            if n == 16:
                near(fd_mean, mean_rise, 0.025, "Independent FD and FEM mean temperature [K]")
                near(fd_center, center_rise, 0.025,
                     "Independent FD and FEM center temperature rise [K]")

                # Compare every exported P1 node with the independent spectral FD field.
                # The six-face Dirichlet boundary is zero rise in both formulations.
                require(len(points) == (n + 1)**3, "Expected a complete Cartesian cube node grid")
                grid_nodes = set()
                nodal_errors = []
                h = length / n
                for point, actual in zip(points, temperature):
                    ijk = tuple(round(coordinate / h) for coordinate in point)
                    require(all(0 <= index <= n for index in ijk),
                            f"NOVA node lies outside the independent FD grid: {point}")
                    for coordinate, index in zip(point, ijk):
                        near(coordinate, index * h, 2e-12,
                             "NOVA node must coincide with an independent FD grid point [m]")
                    require(ijk not in grid_nodes, f"Duplicate independent FD grid node {ijk}")
                    grid_nodes.add(ijk)
                    if any(index in (0, n) for index in ijk):
                        expected = boundary
                    else:
                        expected = boundary + float(fd_rise[ijk[0] - 1, ijk[1] - 1, ijk[2] - 1])
                    nodal_errors.append(actual - expected)
                require(len(grid_nodes) == (n + 1)**3, "Independent FD comparison missed cube nodes")
                maximum_nodal_error = max(abs(error) for error in nodal_errors)
                rms_nodal_error = math.sqrt(math.fsum(error * error for error in nodal_errors) /
                                            len(nodal_errors))
                near(maximum_nodal_error, 0, 2e-8,
                     "Maximum full-field NOVA/SciPy nodal temperature difference [K]")
                near(rms_nodal_error, 0, 2e-8,
                     "RMS full-field NOVA/SciPy nodal temperature difference [K]")
                print(f"full_field_nodes={len(nodal_errors)} max_error_K={maximum_nodal_error:.6g} "
                      f"rms_error_K={rms_nodal_error:.6g}")
            print(f"SciPy FD n={n} mean_rise_K={fd_mean:.12g} center_rise_K={fd_center:.12g} "
                  f"error_upper_K={upper-fd_mean:.8g}")
            previous_mean = fd_mean
            previous_center = fd_center
        require(upper - previous_mean < 0.002, "Finest FD mean exceeds its 0.002 K allowance")


def persistence(runner, fixture, work):
    src, dst = work / "input.nova", work / "saved model.nova"
    vtk = work / "result.vtk"
    base = fixture.read_text(encoding="ascii")
    src.write_bytes(base.replace("\n", "\r\n").encode("ascii"))
    invoke(runner, "save-as", src, dst)
    canonical = dst.read_bytes()
    invoke(runner, "save-as", dst, dst)
    require(dst.read_bytes() == canonical, "Canonical in-place save changed content")
    checked = invoke(runner, "check", dst)
    require("nodes 108\n" in checked and "tetrahedra 288\n" in checked, "Box mesh check")
    first = invoke(runner, "solve", src, "--vtk", vtk)
    output = vtk.read_bytes()
    require(invoke(runner, "solve", dst, "--vtk", vtk) == first and vtk.read_bytes() == output,
            "Save/reload changed numerical reconstruction")
    bad_records = [
        ("NOVA_PROJECT 1", "NOVA_PROJECT 2"), ("project convecting-wall", "project wall"),
        ("tets wall structured_tet", "tets missing structured_tet"),
        ("steady_heat tets homogeneous", "steady_heat tets missing"),
        ("source steady 0", "source other 0"), ("source steady 0", "source steady 0 extra"),
        ("isotropic 2", "isotropic nan"), ("box 0.5", "box 1e9999"),
        ("8 3 2", "-1 3 2"), ("8 3 2", "200000 1 1"), ("8 3 2", "0 3 2"),
        ("8 3 2", "8x 3 2"), ("cg jacobi", "cg unknown"), ("cg jacobi 2000", "cg jacobi 100001"),
        ("xmin temperature 400", "xnone temperature 400"), ("wall xmax", "other xmax"),
        ("convection 12 300", "convection -1 300"), ("temperature 400", "temperature -1"),
        ("end\n", "boundary steady wall xmin temperature 400\nend\n"),
        ("end\n", ""), ("end\n", "end extra\n"), ("end\n", "end\nproject other\n"),
    ]
    for a, b in bad_records:
        require(a in base, f"Missing test mutation {a}")
        src.write_text(base.replace(a, b), encoding="ascii")
        invoke(runner, "save-as", src, dst, success=False)
        require(dst.read_bytes() == canonical, "Invalid input replaced saved project")
    for content in (b"#" * 65537, b"NOVA_PROJECT 1\x00\n", base.replace(
            "box 0.5 0.5 0.5", "box 1 1 1e-50").encode("ascii")):
        src.write_bytes(content)
        invoke(runner, "check", src, success=False)
    src.write_text(base.replace("jacobi 2000", "jacobi 1"), encoding="ascii")
    invoke(runner, "solve", src, "--vtk", vtk, success=False)
    require(vtk.read_bytes() == output, "Nonconverged solve replaced previous VTK")
    src.write_text(base, encoding="ascii")
    as_vtk = work / "input.vtk"
    as_vtk.write_bytes(src.read_bytes())
    invoke(runner, "solve", as_vtk, "--vtk", as_vtk, success=False)
    require(as_vtk.read_bytes() == src.read_bytes(), "Result overwrote its input")
    print("Canonical round-trip, repeatable solve, invalid records and preserved destinations verified.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("suite", choices=("wall", "generated", "cube", "persistence"))
    parser.add_argument("runner", type=Path)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("--scipy", action="store_true", help="Also solve the cube with SciPy finite differences")
    args = parser.parse_args()
    if args.scipy and args.suite != "cube":
        parser.error("--scipy applies to the cube reference")
    with tempfile.TemporaryDirectory() as directory:
        inputs = (args.runner.resolve(), args.fixture.resolve(), Path(directory))
        if args.suite == "cube":
            cube(*inputs, check_scipy=args.scipy)
        else:
            {"wall": wall, "generated": generated, "persistence": persistence}[args.suite](*inputs)
