"""Reproduce the data-driven package study, verify export, conservation and restart."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "examples/package_3d/power_module.nova"


def run(exe, model, folder, *args, error=None):
    p = subprocess.run([str(exe), str(model), str(folder), *map(str, args)],
                       capture_output=True, text=True)
    if error:
        if p.returncode == 0 or error.lower() not in p.stderr.lower():
            raise RuntimeError(f"Expected {error!r}, received {p.stderr!r}")
    elif p.returncode:
        raise RuntimeError(f"{model}: {p.stderr}")
    return p


def frame(path):
    """Read all tuple counts and associations, including multicomponent scalars."""
    lines = iter(path.read_text().splitlines())
    assert next(lines) == "# vtk DataFile Version 3.0"
    next(lines)
    assert next(lines) == "ASCII" and next(lines) == "DATASET UNSTRUCTURED_GRID"
    words = iter(" ".join(lines).split())

    def expect(word):
        val = next(words)
        if val != word:
            raise ValueError(f"Expected {word}, read {val}")

    def vals(count):
        data = [float(next(words)) for _ in range(count)]
        assert all(map(math.isfinite, data))
        return data

    expect("FIELD")
    expect("FieldData")
    for _ in range(int(next(words))):
        next(words)
        a, b = int(next(words)), int(next(words))
        expect("double")
        vals(a*b)
    expect("POINTS")
    n = int(next(words))
    expect("double")
    points = [vals(3) for _ in range(n)]
    expect("CELLS")
    ne = int(next(words))
    assert int(next(words)) == 5*ne
    cells = []
    for _ in range(ne):
        expect("4")
        ids = tuple(int(next(words)) for _ in range(4))
        assert len(set(ids)) == 4 and min(ids) >= 0 and max(ids) < n
        cells.append(ids)
    expect("CELL_TYPES")
    assert int(next(words)) == ne
    for _ in cells:
        expect("10")
    result = {"points": points, "cells": cells}
    count = 0
    for token in words:
        if token in ("POINT_DATA", "CELL_DATA"):
            count = int(next(words))
            assert count == (n if token == "POINT_DATA" else ne)
            continue
        name = next(words)
        assert next(words) in ("double", "int")
        if token == "SCALARS":
            width = int(next(words))
            expect("LOOKUP_TABLE")
            expect("default")
        else:
            width = {"VECTORS": 3, "TENSORS": 9}[token]
        assert name not in result
        result[name] = [vals(width) for _ in range(count)]
    return result


def history(directory):
    with (directory / "history.csv").open(newline="") as f:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(f)]


def load_last(directory):
    collection = ET.parse(directory / "package.pvd").find("Collection")
    times = [float(v.attrib["timestep"]) for v in collection]
    assert times == sorted(set(times))
    return times, frame(directory/collection[-1].attrib["file"])


def verify(exe, work, independent=False):
    full = work/"full"
    run(exe, MODEL, full)
    h = history(full)
    times, f = load_last(full)
    assert times == [r["time_s"] for r in h] and times[-1] == .08
    assert all(t in times for t in [.01, .05, .06])
    assert all(r["error_ratio"] <= 1 for r in h)
    assert sum(r["rejected"] for r in h) > 0
    assert max(abs(r["heat_defect_W"]) for r in h) < 1e-8
    assert max(abs(r["terminal_W"]-r["joule_W"]) for r in h) < 1e-7
    assert max(abs(r["net_current_A"]) for r in h) < 2e-7
    for name in ("temperature", "electric_potential", "displacement", "support_reaction", "electrically_active"):
        assert len(f[name]) == len(f["points"])
    for name in ("material_id", "current_density", "heat_flux", "joule_density", "electrical_conductivity",
                 "young_modulus", "strain", "mean_stress", "von_mises", "principal_stresses"):
        assert len(f[name]) == len(f["cells"])
    assert set(a[0] for a in f["material_id"]) == set(range(1, 7))
    for i, [tag] in enumerate(f["material_id"]):
        if tag in (2, 3, 6):
            assert f["joule_density"][i][0] == 0 and f["current_density"][i] == [0, 0, 0]
        s1, s2, s3 = f["principal_stresses"][i]
        vm = math.sqrt(((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/2)
        assert s1 >= s2-1e-6 and s2 >= s3-1e-6
        assert math.isclose(f["von_mises"][i][0], vm, rel_tol=1e-9, abs_tol=1e-5)
    assert h[-1]["Tmax_K"] > 301 and h[-1]["max_von_mises_Pa"] > 1e6
    assert h[-1]["max_u_m"] > 1e-9
    with (full/"terminals.csv").open() as file:
        currents = list(csv.DictReader(file))
    latest = {int(r["terminal_tag"]): float(r["current_into_domain_A"]) for r in currents if float(r["time_s"]) == .08}
    assert abs(latest[11]) < 1e-12 and abs(latest[12]+latest[13]) < 2e-7
    assert math.isclose(.75*latest[13], h[-1]["joule_W"], abs_tol=1e-7)

    cut = times.index(.05)+1
    checkpoint = full/f"step_{cut:06}.chk"
    restarted = work/"restart"
    run(exe, MODEL, restarted, "--restart", checkpoint)
    suffix = history(restarted)
    assert suffix == h[cut:], "Restart changed the multistep/controller trajectory"
    assert load_last(restarted)[1] == f
    run(exe, MODEL, full, error="already exists")
    damaged = work/"damaged.chk"
    damaged.write_bytes(checkpoint.read_bytes()[:-25])
    run(exe, MODEL, work/"broken", "--restart", damaged, error="checkpoint")
    altered = work/"altered.nova"
    # Absolute mesh path lets the mutated model live in the temporary directory.
    altered.write_text(MODEL.read_text().replace("mesh package.msh", f"mesh {MODEL.with_name('package.msh').as_posix()}")
                       .replace("adapt 0.015", "adapt 0.02"))
    run(exe, altered, work/"mismatch", "--restart", checkpoint, error="identity")

    result = json.loads((full/"summary.json").read_text())
    if independent:
        from reference_package_3d import compare, mechanical_and_transient
        result["independent_electrical"] = compare(f, latest[13])
        result["independent_mechanics_and_heat"] = mechanical_and_transient(f)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--independent", action="store_true", help="also solve DC equilibrium independently with SciPy")
    args = parser.parse_args()
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        print(json.dumps(verify(args.executable.resolve(), args.output, args.independent), indent=2))
    else:
        with tempfile.TemporaryDirectory() as d:
            print(json.dumps(verify(args.executable.resolve(), Path(d), args.independent), indent=2))
