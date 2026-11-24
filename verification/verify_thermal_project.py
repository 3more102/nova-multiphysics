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
    parser.add_argument("suite", choices=("wall", "generated", "persistence"))
    parser.add_argument("runner", type=Path)
    parser.add_argument("fixture", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as directory:
        {"wall": wall, "generated": generated, "persistence": persistence}[args.suite](
            args.runner.resolve(), args.fixture.resolve(), Path(directory))
