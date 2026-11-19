"""End-to-end steady conduction verification using only the Python standard library."""

import argparse
import itertools
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import tempfile
import time


def write_wall(path, nx, ny, nz):
    """Conforming cube with shared interface vertices, sparse IDs and mixed windings."""
    def node(i, j, k):
        return 10 + 7 * (i + (nx + 1) * (j + (ny + 1) * k))

    points = [(node(i, j, k), i / nx, j / ny, k / nz)
              for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    elements = []
    for i, tag in ((0, 101), (nx, 102)):
        for j, k in itertools.product(range(ny), range(nz)):
            a, b = node(i, j, k), node(i, j + 1, k)
            c, d = node(i, j + 1, k + 1), node(i, j, k + 1)
            elements.extend(((2, tag, (a, b, c)), (2, tag, (a, c, d))))
    for k, j, i in itertools.product(range(nz), range(ny), range(nx)):
        tag = 11 if 2 * i < nx else 22
        for axes in itertools.permutations(range(3)):
            vertex = [i, j, k]
            tetra = [node(*vertex)]
            for axis in axes:
                vertex[axis] += 1
                tetra.append(node(*vertex))
            elements.append((4, tag, tetra))
    with path.open("w", encoding="ascii") as out:
        out.write(f"$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n{len(points)}\n")
        for ident, x, y, z in reversed(points):
            out.write(f"{ident} {x:.17g} {y:.17g} {z:.17g}\n")
        out.write(f"$EndNodes\n$Elements\n{len(elements)}\n")
        for index, (kind, tag, nodes) in enumerate(elements, 1):
            out.write(f"{index} {kind} 2 {tag} 1 " + " ".join(map(str, nodes)) + "\n")
        out.write("$EndElements\n")


def read_fields(path, require_heat_generation=False):
    """Read exactly the numerical fields that a downstream VTK consumer receives."""
    lines = path.read_text(encoding="ascii").splitlines()
    if lines[0] != "# vtk DataFile Version 3.0" or lines[2] != "ASCII":
        raise ValueError("Invalid VTK header")
    tokens = iter(" ".join(lines[3:]).split())

    def expect(words):
        for word in words.split():
            if next(tokens) != word:
                raise ValueError(f"Expected VTK token {word}")

    def real():
        value = float(next(tokens))
        if not math.isfinite(value):
            raise ValueError("Nonfinite VTK data")
        return value

    expect("DATASET UNSTRUCTURED_GRID POINTS")
    n = int(next(tokens))
    expect("double")
    points = [tuple(real() for _ in range(3)) for _ in range(n)]
    expect("CELLS")
    m = int(next(tokens))
    expect(str(5 * m))
    cells = []
    for _ in range(m):
        expect("4")
        cell = tuple(int(next(tokens)) for _ in range(4))
        if min(cell) < 0 or max(cell) >= n or len(set(cell)) != 4:
            raise ValueError("Invalid VTK tetrahedral connectivity")
        cells.append(cell)
    expect(f"CELL_TYPES {m}")
    for _ in cells:
        expect("10")
    expect(f"POINT_DATA {n} SCALARS temperature double 1 LOOKUP_TABLE default")
    temperature = [real() for _ in points]
    expect(f"CELL_DATA {m} SCALARS material_id int 1 LOOKUP_TABLE default")
    tags = [int(next(tokens)) for _ in cells]
    token = next(tokens)
    sources = None
    if token == "SCALARS":
        expect("heat_generation double 1 LOOKUP_TABLE default")
        sources = [real() for _ in cells]
        token = next(tokens)
    if token != "VECTORS":
        raise ValueError("Missing VTK heat flux vectors")
    expect("heat_flux double")
    flux = [tuple(real() for _ in range(3)) for _ in cells]
    if next(tokens, None) is not None:
        raise ValueError("Unexpected trailing VTK data")
    for cell, tag in zip(cells, tags):
        xs = [points[i][0] for i in cell]
        if not ((tag == 11 and max(xs) <= 0.5) or (tag == 22 and min(xs) >= 0.5)):
            raise ValueError("Exported volume ID does not match its geometric layer")
    if require_heat_generation:
        if sources is None:
            raise ValueError("Expected cell heat_generation field")
        return points, cells, temperature, flux, sources
    return points, cells, temperature, flux


MATERIALS = ["--material", "22", "6", "--material", "11", "2"]
HOT = ["--temperature", "101", "400"]
ROBIN = ["--convection", "102", "12", "300"]


def measure(label, actual, expected, atol):
    error = abs(actual - expected)
    relative = f"{error / abs(expected):.4e}" if expected else "n/a"
    passed = math.isfinite(actual) and error <= atol
    print(f"{label}: expected={expected:.12g} computed={actual:.17g} "
          f"abs_error={error:.4e} rel_error={relative} limit={atol:.2e} "
          f"{'PASS' if passed else 'FAIL'}")
    if not passed:
        raise ValueError(label)


def run(exe, mesh, vtk, args, check=True):
    return subprocess.run([exe, mesh, *args, "--vtk", vtk],
                          capture_output=True, text=True, check=check)


def read_balance(stdout):
    """Parse independently printed global FEM energy and reaction diagnostics."""
    keys = {
        "FEM generated heat [W]": "source",
        "FEM Neumann outward heat [W]": "neumann",
        "FEM Robin outward heat [W]": "robin",
        "FEM constrained-node outward reaction [W]": "constraint",
        "FEM signed energy balance defect [W]": "defect",
        "FEM maximum free-equation residual [W]": "free",
    }
    fields = {}
    for line in stdout.splitlines():
        for prefix, name in keys.items():
            if line.startswith(prefix + ": "):
                if name in fields:
                    raise ValueError("Repeated FEM energy-balance field: " + prefix)
                fields[name] = float(line.split(": ", 1)[1])
    if set(fields) != set(keys.values()) or not all(map(math.isfinite, fields.values())):
        raise ValueError("FEM energy-balance report is missing fields or non-finite")
    return fields


def check_wall(exe, mesh, vtk):
    execution = run(exe, mesh, vtk, MATERIALS + HOT + ROBIN
                    + ["--flux-report", "--balance-report"])
    energy = read_balance(execution.stdout)
    measure("FEM wall generation [W]", energy["source"], 0, 2e-8)
    measure("FEM wall Robin heat [W]", energy["robin"], 240, 2e-7)
    measure("FEM wall Dirichlet reaction [W]", energy["constraint"], -240, 2e-7)
    measure("FEM wall energy-balance defect [W]", energy["defect"], 0, 2e-7)
    measure("FEM wall maximum free residual [W]", energy["free"], 0, 2e-7)
    reported = {}
    for line in execution.stdout.splitlines():
        if line.startswith("  physical tag "):
            label, rest = line.strip().split(": faces=", 1)
            count, rest = rest.split(" area_m2=", 1)
            area, rate = rest.split(" outward_W=", 1)
            tag = int(label.removeprefix("physical tag "))
            if tag in reported:
                raise ValueError("Duplicate boundary-report tag")
            reported[tag] = (int(count), float(area), float(rate))
    if set(reported) != {0, 101, 102}:
        raise ValueError(f"Missing or unexpected boundary tags: {reported}")
    if reported[101][0] != 2 or reported[102][0] != 2:
        raise ValueError("Unexpected number of tagged wall boundary faces")
    measure("tagged hot face area [m2]", reported[101][1], 1, 1e-12)
    measure("tagged cold face area [m2]", reported[102][1], 1, 1e-12)
    measure("tagged hot outward rate [W]", reported[101][2], -240, 2e-7)
    measure("tagged cold outward rate [W]", reported[102][2], 240, 2e-7)
    measure("untagged insulated sides [W]", reported[0][2], 0, 2e-7)
    net_line = next((line for line in execution.stdout.splitlines()
                     if line.startswith("Net outward conductive heat rate [W]: ")), None)
    if net_line is None:
        raise ValueError("No net conductive heat-rate report")
    measure("reported net outward rate [W]", float(net_line.split(": ", 1)[1]), 0, 2e-7)
    points, cells, temperatures, fluxes = read_fields(vtk)
    # q = (400 - 300)/(0.5/2 + 0.5/6 + 1/12); both layers have unit area.
    measure("maximum nodal temperature error [K]", max(
        abs(t - (400 - 120 * min(x, 0.5) - 40 * max(x - 0.5, 0.0)))
        for (x, _, _), t in zip(points, temperatures)), 0, 2e-8)
    for x, target in ((0.5, 340), (1.0, 320)):
        surface = [t for p, t in zip(points, temperatures) if p[0] == x]
        measure(f"worst temperature at x={x} [K]", max(surface, key=lambda t: abs(t-target)),
                target, 2e-8)
    measure("worst axial heat flux [W/m2]", max((q[0] for q in fluxes),
            key=lambda v: abs(v - 240)), 240, 2e-7)
    measure("maximum transverse flux [W/m2]", max(math.hypot(y, z) for _, y, z in fluxes),
            0, 2e-7)
    rates = [0.0, 0.0]
    for cell, (qx, _, _) in zip(cells, fluxes):
        for x in (0, 1):
            face = [points[node] for node in cell if points[node][0] == x]
            if len(face) == 3:
                a, b, c = face
                area = 0.5 * abs((b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]))
                rates[x] += (2*x - 1) * qx * area
    measure("hot-face outward rate [W]", rates[0], -240, 2e-7)
    measure("cold-face outward rate [W]", rates[1], 240, 2e-7)
    measure("net outward rate [W]", sum(rates), 0, 2e-7)


def quadrature():
    # Duffy map of [0,1]^3 to the reference tetrahedron. For a degree-four
    # polynomial the transformed degrees are at most 6, 5, 4, respectively.
    g4 = [(0.06943184420297371, 0.1739274225687269), (0.33000947820757187, 0.3260725774312731),
          (0.6699905217924281, 0.3260725774312731), (0.9305681557970262, 0.1739274225687269)]
    a = math.sqrt(3/5) / 2
    g3 = [(0.5 - a, 5/18), (0.5, 4/9), (0.5 + a, 5/18)]
    for (u, wu), (v, wv), (w, ww) in itertools.product(g4, g3, g3):
        r, s, t = u, (1-u)*v, (1-u)*(1-v)*w
        yield (1-r-s-t, r, s, t), wu*wv*ww*(1-u)**2*(1-v)


def field_errors(fields, n):
    points, cells, temperatures, fluxes = fields
    te = qe = 0.0
    samples = list(quadrature())
    # Each grid cell decomposes into six equal-volume tetrahedra, |det J| = 1/n^3.
    for cell, (qx, qy, qz) in zip(cells, fluxes):
        xs = [points[i][0] for i in cell]
        ts = [temperatures[i] for i in cell]
        for shape, w in samples:
            x = sum(b * v for b, v in zip(shape, xs))
            t = sum(b * v for b, v in zip(shape, ts))
            primitive = 120*x + 60*x*x
            exact = 400 - primitive/2 if max(xs) <= 0.5 else 362.5 - (primitive - 75)/6
            te += w * (t - exact)**2
            qe += w * ((qx - (120 + 120*x))**2 + qy*qy + qz*qz)
    return math.sqrt(te / n**3), math.sqrt(qe / n**3)


def check_refinement(exe, directory):
    previous = None
    for n in (2, 4, 8):
        mesh, vtk = directory / f"source-{n}.msh", directory / f"source-{n}.vtk"
        write_wall(mesh, n, n, n)
        run(exe, mesh, vtk, MATERIALS + HOT + ["--source", "120", "--flux", "102", "240"])
        errors = field_errors(read_fields(vtk), n)
        print(f"n={n} tetrahedra={6*n**3} temperature_L2={errors[0]:.10g} flux_L2={errors[1]:.10g}")
        if previous:
            orders = [math.log2(old/new) for old, new in zip(previous, errors)]
            # Allow the pre-asymptotic coarse mesh to fall below the nominal P1 orders.
            print(f"temperature_order={orders[0]:.8f} flux_order={orders[1]:.8f}")
            if not (orders[0] >= 1.7 and orders[1] >= 0.85):
                raise ValueError("Insufficient refinement rate")
        previous = errors
    if previous[0] > 0.08 or previous[1] > 5.0:
        raise ValueError("Fine-grid temperature or flux error exceeds its absolute limit")


def check_tagged_sources(exe, directory):
    """Two-layer wall with distinct internal heating; independent 1D reference."""
    previous = None
    for n in (2, 4, 8):
        mesh, vtk = directory / f"tagged-heat-{n}.msh", directory / f"tagged-heat-{n}.vtk"
        write_wall(mesh, n, n, n)
        execution = run(exe, mesh, vtk, MATERIALS + HOT + ROBIN
            + ["--heat-source", "11", "120", "--heat-source", "22", "360",
               "--balance-report"])
        energy = read_balance(execution.stdout)
        measure(f"FEM regional generation n={n} [W]", energy["source"], 240, 2e-7)
        measure(f"FEM regional balance defect n={n} [W]", energy["defect"], 0, 2e-7)
        measure(f"FEM regional max free residual n={n} [W]", energy["free"], 0, 2e-7)
        measure(f"FEM regional summed outward n={n} [W]",
                energy["neumann"] + energy["robin"] + energy["constraint"],
                240, 2e-7)
        points, cells, temperatures, fluxes, sources = read_fields(
            vtk, require_heat_generation=True)
        for cell, source in zip(cells, sources):
            xs = [points[node][0] for node in cell]
            expected = 120 if max(xs) <= 0.5 else 360
            if source != expected:
                raise ValueError("Exported heat generation does not match physical volume")
        te = qe = 0.0
        for cell, (qx, qy, qz) in zip(cells, fluxes):
            xs = [points[i][0] for i in cell]
            ts = [temperatures[i] for i in cell]
            layer_one = max(xs) <= 0.5
            for shape, weight in quadrature():
                x = sum(b * value for b, value in zip(shape, xs))
                t = sum(b * value for b, value in zip(shape, ts))
                if layer_one:
                    exact_t = 400 - 72*x - 30*x*x
                    exact_q = 144 + 120*x
                else:
                    distance = x - 0.5
                    exact_t = 356.5 - 34*distance - 30*distance*distance
                    exact_q = 204 + 360*distance
                te += weight * (t - exact_t)**2
                qe += weight * ((qx - exact_q)**2 + qy*qy + qz*qz)
        errors = (math.sqrt(te/n**3), math.sqrt(qe/n**3))
        print(f"tagged sources n={n} T_L2={errors[0]:.10g} flux_L2={errors[1]:.10g}")
        if previous is not None:
            orders = [math.log2(old / new) for old, new in zip(previous, errors)]
            print(f"tagged temperature_order={orders[0]:.8f} flux_order={orders[1]:.8f}")
            if not (orders[0] >= 1.65 and orders[1] >= 0.80):
                raise ValueError("Insufficient convergence for tagged internal heat generation")
        previous = errors
        if n == 2:
            # A uniform spatial source is additive to each physical-region source.
            run(exe, mesh, vtk, MATERIALS + HOT + ROBIN +
                ["--source", "5", "--heat-source", "11", "120",
                 "--heat-source", "22", "360"])
            _, _, _, _, effective = read_fields(vtk, require_heat_generation=True)
            for cell, source in zip(cells, effective):
                expected = 125 if max(points[node][0] for node in cell) <= 0.5 else 365
                if source != expected:
                    raise ValueError("VTK heat generation omitted uniform source contribution")
        if n == 8:
            cold = [t for (x, _, _), t in zip(points, temperatures) if x == 1]
            interface = [t for (x, _, _), t in zip(points, temperatures) if x == 0.5]
            measure("tagged heat interface [K]", max(interface, key=lambda t: abs(t-356.5)),
                    356.5, 0.5)
            measure("tagged heat convection surface [K]", max(cold, key=lambda t: abs(t-332)),
                    332, 0.5)
            # Integrated sources: 120*0.5 + 360*0.5 = 240 W on unit area.
            # The 332 K cold surface gives 12*(332-300)=384 W convection;
            # hence the hot boundary supplies the remaining 144 W.
            measure("tagged heat convection rate [W]", 12*(sum(cold)/len(cold)-300),
                    384, 6.0)
    if previous[0] > 0.4 or previous[1] > 18.0:
        raise ValueError("Tagged internal-heating fine-grid errors exceed limits")


def check_invalid_cases(exe, mesh, vtk):
    # Reverse the material assignments without changing connectivity. The total
    # resistance stays fixed while the interface temperature moves by 40 K.
    run(exe, mesh, vtk, ["--material", "11", "6", "--material", "22", "2"] + HOT + ROBIN)
    points, _, temperature, _ = read_fields(vtk)
    interface = [t for p, t in zip(points, temperature) if p[0] == 0.5]
    measure("swapped material interface [K]", max(interface, key=lambda t: abs(t-380)), 380, 2e-8)
    original = vtk.read_bytes()
    cases = [
        (["--material", "11", "2"] + HOT, "Missing conductivity"),
        (MATERIALS + ["--material", "7", "1"], "No volume"),
        (MATERIALS + HOT + ["--flux", "101", "3"], "Duplicate surface"),
        (MATERIALS + ["--temperature", "777", "0"], "No exterior face"),
        (["--material", "11", "0"], "positive"),
        (MATERIALS + ["--material", "22", "3"], "Duplicate material"),
        (MATERIALS + ["--source", "nan"], "finite number"),
        (MATERIALS + HOT + ["--balance-report", "--balance-report"], "Duplicate --balance-report"),
        (MATERIALS + HOT + ["--heat-source", "777", "1"], "No volume"),
        (MATERIALS + HOT + ["--heat-source", "11", "1", "--heat-source", "11", "2"],
         "Duplicate heat source"),
        (MATERIALS + HOT + ["--heat-source", "22", "inf"], "finite number"),
        (MATERIALS + ["--heat-source", "11", "10"], "Incompatible pure-Neumann"),
        (MATERIALS + ["--source", "1"], "Incompatible pure-Neumann"),
        (MATERIALS + ["--convection", "102", "-1", "300"], "nonnegative"),
    ]
    for args, message in cases:
        result = run(exe, mesh, vtk, args, check=False)
        if result.returncode == 0 or message not in result.stderr or vtk.read_bytes() != original:
            raise ValueError(f"Invalid case mishandled: {args}: {result.stderr}")
    before = mesh.read_bytes()
    result = run(exe, mesh, mesh, MATERIALS + HOT, check=False)
    if result.returncode == 0 or "overwrite" not in result.stderr or mesh.read_bytes() != before:
        raise ValueError("Mesh/output alias mishandled")
    print(f"{len(cases)+1} invalid CLI cases rejected with input and prior output preserved")


def benchmark(exe, directory):
    for n in (8, 16, 24):
        mesh, vtk = directory / "bench.msh", directory / "bench.vtk"
        write_wall(mesh, n, n, n)
        elapsed = []
        for _ in range(3):
            start = time.perf_counter()
            result = run(exe, mesh, vtk, MATERIALS + HOT + ROBIN)
            elapsed.append(time.perf_counter() - start)
        iterations = next(line for line in result.stdout.splitlines() if line.startswith("CG iterations"))
        print(f"n={n} DOFs={(n+1)**3} tetrahedra={6*n**3} {iterations} "
              f"median_wall_s={statistics.median(elapsed):.6f} vtk_bytes={vtk.stat().st_size}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exe", type=Path)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("--benchmark", action="store_true")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as work:
        root = Path(work)
        output = root / "wall.vtk"
        check_wall(args.exe.resolve(), args.fixture, output)
        private_fixture = root / "input.msh"
        shutil.copyfile(args.fixture, private_fixture)
        check_invalid_cases(args.exe.resolve(), private_fixture, output)
        check_refinement(args.exe.resolve(), root)
        check_tagged_sources(args.exe.resolve(), root)
        if args.benchmark:
            benchmark(args.exe.resolve(), root)
