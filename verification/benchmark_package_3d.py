"""Run each deterministic workload in a fresh process and emit one JSON record per run."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


def measure(exe, n, mode, timeout, memory_mib):
    def limit():
        import resource
        bound = memory_mib * 1024**2
        resource.setrlimit(resource.RLIMIT_AS, (bound, bound))

    start = time.monotonic()
    rss = None
    with tempfile.TemporaryFile(mode="w+t") as output:
        with subprocess.Popen([str(exe), str(n), mode], stdout=output, stderr=subprocess.STDOUT,
                              preexec_fn=limit if memory_mib and os.name == "posix" else None) as proc:
            if hasattr(os, "wait4"):
                while True:
                    pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
                    if pid:
                        proc.returncode = os.waitstatus_to_exitcode(status)
                        rss = usage.ru_maxrss / (1024**2 if sys.platform == "darwin" else 1024)
                        break
                    if time.monotonic()-start > timeout:
                        proc.kill()
                    time.sleep(.01)
            else:
                try:
                    proc.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
            output.seek(0)
            return dict(subdivisions=n, mode=mode, exit=proc.returncode,
                        wall_s=round(time.monotonic()-start, 5), peak_rss_MiB=rss,
                        result=output.read().strip())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--sizes", type=int, nargs="+", default=[12, 27, 46, 58])
    parser.add_argument("--mode", choices=["assemble", "solve"], default="solve")
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--memory-mib", type=int, default=0)
    args = parser.parse_args()
    failures = 0
    for n in args.sizes:
        row = measure(args.executable.resolve(), n, args.mode, args.timeout, args.memory_mib)
        print(json.dumps(row), flush=True)
        failures += row["exit"] != 0
    sys.exit(bool(failures))
