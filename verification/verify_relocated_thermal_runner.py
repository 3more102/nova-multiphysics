"""Install a built shared SDK, move its prefix and consume it from a clean directory."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def one(prefix, name):
    found = list(prefix.rglob(name))
    if len(found) != 1:
        raise RuntimeError(f"Expected exactly one {name!r} under the relocated prefix, found {found}")
    return found[0]


def verify(build, cmake, ctest):
    repository = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        initial = root / "original prefix"
        moved = root / "moved SDK"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        env.pop("LD_PRELOAD", None)
        env.pop("NOVA_DIR", None)

        def run(args, **kwargs):
            return subprocess.run([str(x) for x in args], cwd=root, env=env, check=True, **kwargs)

        run([cmake, "--install", build, "--config", "Release", "--prefix", initial])
        initial.rename(moved)
        runner = one(moved, "nova-run.exe" if os.name == "nt" else "nova-run")
        project = one(moved, "convecting_wall.nova")
        if sys.platform.startswith("linux"):
            core = one(moved, "libnova_core.so")
            deps = run(["ldd", runner], text=True, capture_output=True).stdout
            resolved = re.search(r"libnova_core\.so\s*=>\s*(.+?)\s+\(0x[0-9a-fA-F]+\)", deps)
            if not resolved or not Path(resolved[1]).samefile(core):
                raise RuntimeError(f"Runner did not load the relocated NOVA library:\n{deps}")
        elif os.name == "nt":
            dlls = list(runner.parent.glob("*nova_core.dll"))
            if len(dlls) != 1:
                raise RuntimeError(f"Expected the NOVA runtime DLL alongside the solver, found {dlls}")

        # A stand-alone runner must locate its own library without modifying PATH.
        run([sys.executable, repository / "verification/verify_thermal_project.py", "wall", runner, project])
        if os.name == "nt":
            # External consumer executables are outside the SDK bin directory.
            env["PATH"] = str(runner.parent) + os.pathsep + env.get("PATH", "")
        installed_config = one(moved, "NOVAConfig.cmake").parent
        consumer = root / "consumer build"
        run([cmake, "-S", repository / "tests/package_consumer", "-B", consumer,
             "-DCMAKE_BUILD_TYPE=Release", f"-DNOVA_DIR={installed_config}"])
        run([cmake, "--build", consumer, "--config", "Release", "--parallel", "3"])
        run([ctest, "--test-dir", consumer, "-C", "Release", "--output-on-failure"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--ctest", default="ctest")
    args = parser.parse_args()
    verify(args.build.resolve(), args.cmake, args.ctest)
