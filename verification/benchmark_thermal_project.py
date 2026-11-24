"""Linux single-process wall benchmark; measures one run per requested grid size."""
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from time import perf_counter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runner", type=Path)
    parser.add_argument("fixture", type=Path)
    parser.add_argument("sizes", type=int, nargs="*", default=[8, 16, 32])
    args = parser.parse_args()
    if not sys.platform.startswith("linux"):
        parser.error("Peak RSS in this benchmark uses Linux wait4 KiB units.")
    fixture = args.fixture.read_text(encoding="ascii")
    if "structured_tet 8 3 2" not in fixture:
        parser.error("Expected the distributed convecting_wall.nova reference fixture.")
    csvout = csv.writer(sys.stdout, lineterminator="\n")
    csvout.writerow(("divisions", "nodes", "tetrahedra", "wall_s", "peak_rss_KiB",
                     "cg_iterations", "cg_residual", "cold_error_K", "convection_error_W", "balance_W"))
    with tempfile.TemporaryDirectory() as directory:
        project = Path(directory) / "wall.nova"
        for n in args.sizes:
            project.write_text(fixture.replace("structured_tet 8 3 2", f"structured_tet {n} {n} {n}"),
                               encoding="ascii")
            with tempfile.TemporaryFile(mode="w+") as out:
                start = perf_counter()
                with subprocess.Popen([str(args.runner.resolve()), "solve", str(project)], stdout=out) as child:
                    _, status, usage = os.wait4(child.pid, 0)
                    child.returncode = os.waitstatus_to_exitcode(status)
                    if child.returncode:
                        raise RuntimeError(f"Solver exited with status {child.returncode}")
                elapsed = perf_counter() - start
                out.seek(0)
                data = dict(line.split() for line in out if line.strip())
            cold_error = abs(float(data["min_K"]) - 325.0)
            heat_error = abs(float(data["convection_out_W"]) - 75.0)
            balance = float(data["balance_defect_W"])
            if not all(math.isfinite(x) and abs(x) <= 2e-8 for x in (cold_error, heat_error, balance)):
                raise AssertionError(f"Wall accuracy threshold exceeded on {n}^3 grid")
            csvout.writerow((n, data["nodes"], data["tetrahedra"], f"{elapsed:.6f}", usage.ru_maxrss,
                             data["iterations"], data["residual"], cold_error, heat_error, balance))


if __name__ == "__main__":
    main()
