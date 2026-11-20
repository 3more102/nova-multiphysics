#!/usr/bin/env python3
"""Independent end-to-end checks for tagged transient tetrahedral quality gates."""
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def run(cli, mesh, *args):
    return subprocess.run([str(cli), str(mesh), *map(str, args)],
                          capture_output=True, text=True, check=False)


def metric(output, name):
    match = re.search(r"^" + re.escape(name) + r":\s*([-+0-9.eE]+)\s*$",
                      output, re.MULTILINE)
    if not match:
        raise AssertionError(f"Missing {name} in quality report: {output}")
    return float(match.group(1))


def verify(cli, mesh):
    inspected = run(cli, mesh, "--check-mesh")
    if inspected.returncode:
        raise AssertionError(f"Mesh inspection failed: {inspected.stderr}")
    expected = {
        "Total tetrahedral volume [m^3]": 1.0 / 6.0,
        "Minimum mean ratio": 12.0 * (0.5 ** (2.0 / 3.0)) / 9.0,
        "Minimum dihedral angle [degrees]":
            math.degrees(math.acos(1.0 / math.sqrt(3.0))),
    }
    for name, reference in expected.items():
        actual = metric(inspected.stdout, name)
        if not math.isclose(actual, reference, rel_tol=1e-12, abs_tol=1e-12):
            raise AssertionError(f"{name}: expected {reference}, got {actual}")

    accepted = run(cli, mesh, "--min-mean-ratio", "0.8", "--check-mesh")
    if accepted.returncode:
        raise AssertionError(f"Valid inspection threshold rejected: {accepted.stderr}")
    rejected = run(cli, mesh, "--check-mesh", "--min-mean-ratio", "0.9")
    if rejected.returncode == 0 or "quality gate failed" not in rejected.stderr:
        raise AssertionError("Strict inspection threshold was not enforced")

    thermal = ("--material", "solid", "1",
               "--capacity", "solid", "2",
               "--initial", "300", "--dt", "0.1", "--steps", "1",
               "--temperature", "fixed-temperature", "300")
    for options in (("--quality-report", "--min-mean-ratio", "0.8"),
                    ("--min-mean-ratio", "0.8", "--quality-report"),
                    ("--min-mean-ratio", "0.8")):
        passed = run(cli, mesh, *thermal, *options)
        if passed.returncode or "NOVA Gmsh transient 3D heat study" not in passed.stdout:
            raise AssertionError(f"Valid quality-gated thermal solve failed: {passed.stderr}")
        if "Minimum mean ratio" not in passed.stdout:
            raise AssertionError("Quality-gated solve omitted the quality report")

    with tempfile.TemporaryDirectory(prefix="nova-gmsh-quality-") as folder:
        # Both output modes must reject a bad mesh before publishing anything.
        for option, destination in (
            ("--vtk", Path(folder) / "should_not_exist.vtk"),
            ("--vtk-series", Path(folder) / "should_not_exist.pvd"),
        ):
            blocked = run(cli, mesh, *thermal, "--min-mean-ratio", "0.9",
                          option, destination)
            if (blocked.returncode == 0 or
                    "quality gate failed" not in blocked.stderr or
                    "final_time" in blocked.stdout or destination.exists()):
                raise AssertionError(
                    f"Failed pre-solve quality gate allowed {option} output")
        if list(Path(folder).glob("should_not_exist_*.vtk")):
            raise AssertionError("Failed quality gate left orphan ParaView frames")

    for invalid in ("0", "-0.1", "1.01", "nan", "inf", "not-a-number"):
        bad = run(cli, mesh, "--check-mesh", "--min-mean-ratio", invalid)
        if bad.returncode == 0:
            raise AssertionError(f"Invalid threshold {invalid!r} accepted")

    for options in (("--check-mesh", "--material", "solid", "1"),
                    ("--material", "solid", "1", "--check-mesh"),
                    ("--check-mesh", "--initial", "300"),
                    ("--check-mesh", "--vtk-series", "foo.pvd"),
                    ("--check-mesh", "--scheme", "cn"),
                    ("--check-mesh", "--check-mesh"),
                    ("--quality-report", "--quality-report"),
                    ("--min-mean-ratio", "0.8", "--min-mean-ratio", "0.8")):
        if run(cli, mesh, *options).returncode == 0:
            raise AssertionError(f"Conflicting/duplicate options accepted: {options}")

    print(f"PASS: transient Gmsh mesh quality gates and analytical checks ({mesh.name})")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: verify_gmsh_transient_quality_3d.py cli mesh")
    verify(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
