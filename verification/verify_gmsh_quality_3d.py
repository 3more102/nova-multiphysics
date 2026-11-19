#!/usr/bin/env python3
"""Independent analytical end-to-end check for Gmsh tetrahedral quality gates."""
import math
import re
import subprocess
import sys


def execute(cli, mesh, *args):
    return subprocess.run([cli, mesh, *args], capture_output=True, text=True, check=False)


def value(text, label):
    match = re.search(r"^" + re.escape(label) + r":\s*([-+0-9.eE]+)\s*$", text, re.M)
    if not match:
        raise AssertionError(f"Missing field {label}: {text}")
    return float(match.group(1))


def check():
    if len(sys.argv) != 3:
        raise SystemExit("usage: verify_gmsh_quality_3d.py cli mesh")
    cli, mesh = sys.argv[1:]
    result = execute(cli, mesh, "--check-mesh")
    if result.returncode != 0:
        raise AssertionError(f"Inspection failed: {result.stderr}")
    # Unit right tetrahedron with vertices at the origin and unit axes.
    ratio = 12.0 * (0.5 ** (2.0 / 3.0)) / 9.0
    angle = math.degrees(math.acos(1.0 / math.sqrt(3.0)))
    checks = (
        ("Total tetrahedral volume [m^3]", 1.0 / 6.0),
        ("Minimum mean ratio", ratio),
        ("Minimum dihedral angle [degrees]", angle),
    )
    for label, expected in checks:
        found = value(result.stdout, label)
        if not math.isclose(found, expected, abs_tol=1e-12, rel_tol=1e-12):
            raise AssertionError(f"{label}: expected {expected}, received {found}")
    accepted = execute(cli, mesh, "--check-mesh", "--min-mean-ratio", "0.8")
    if accepted.returncode != 0:
        raise AssertionError(f"Valid quality gate failed: {accepted.stderr}")
    rejected = execute(cli, mesh, "--check-mesh", "--min-mean-ratio", "0.9")
    if rejected.returncode == 0 or "quality gate failed" not in rejected.stderr:
        raise AssertionError("Low-quality mesh was not rejected by strict threshold")
    for invalid in ("0", "-0.1", "1.01", "nan"):
        bad = execute(cli, mesh, "--check-mesh", "--min-mean-ratio", invalid)
        if bad.returncode == 0:
            raise AssertionError(f"Invalid quality gate {invalid} accepted")
    conflicting = execute(cli, mesh, "--check-mesh", "--material", "solid", "1")
    if conflicting.returncode == 0:
        raise AssertionError("Incompatible quality-only and simulation flags accepted")
    print("Gmsh quality report, analytical values, gate and invalid-argument checks passed")


if __name__ == "__main__":
    check()
