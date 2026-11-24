"""Time the matched SDK workload and the complete multi-field Gmsh CLI."""
import argparse
import math
from pathlib import Path
from statistics import median
import subprocess
import tempfile

from verify_tagged_conduction_3d import write_wall


def run(args, timer):
    # Optional GNU time executable (or compatible wrapper) reports elapsed
    # seconds and peak RSS in KiB, once for each fresh solver process.
    if timer:
        args = [timer, "-f", "%e %M", *args]
    p = subprocess.run(list(map(str,args)),capture_output=True,text=True,check=True)
    rss = int(p.stderr.splitlines()[-1].split()[1])/1024 if timer else None
    return p.stdout, rss


def benchmark(args, directory):
    from time import perf_counter
    for n in args.sizes:
        mesh = directory/f"mesh{n}.msh"
        write_wall(mesh,n,n,n)
        modes = {"prepared": args.prepared}
        if args.reference:
            modes["original"] = args.reference
        history = {mode: [] for mode in [*modes,"cli"]}
        for rep in range(args.repeat):
            matched = {}
            for mode,exe in (list(modes.items()) if rep%2 == 0 else reversed(list(modes.items()))):
                report, rss = run([exe,mesh,directory/f"{mode}.vtk"],args.time_exe)
                words = report.split()
                metrics = {key: float(value) for key,value in zip(words[::2],words[1::2])}
                matched[mode] = metrics
                history[mode].append(dict(metrics, **({"rss_MiB":rss} if rss is not None else {})))
            if args.reference:
                for key in ("energy","displacement_norm","stress_norm"):
                    if not math.isclose(matched["original"][key],matched["prepared"][key],rel_tol=2e-8,abs_tol=1e-12):
                        raise ValueError(f"Unmatched results at n={n}, {key}: {matched}")
            with tempfile.TemporaryDirectory(dir=directory) as frames:
                pvd = Path(frames)/"heat.pvd"
                start = perf_counter()
                _, rss = run([args.cli,mesh,"--material",11,2,"--material",22,2,
                              "--capacity",11,2,"--capacity",22,2,"--initial",310,"--source",-4,
                              "--dt",.125,"--steps",5,"--scheme","cn","--clamp",101,
                              "--elastic",11,1000,.2,.001,300,"--elastic",22,2000,.2,.001,300,
                              "--vtk-series",pvd],args.time_exe)
                elapsed = perf_counter()-start
                if not pvd.is_file() or len(list(Path(frames).glob("*.vtk"))) != 6:
                    raise ValueError("Incomplete VTK series")
                history["cli"].append({"total_s":elapsed,**({"rss_MiB":rss} if rss is not None else {})})
        for mode,runs in history.items():
            stats = " ".join(f"{key}={median(s[key] for s in runs):.8g}" for key in runs[0])
            print(f"n={n} nodes={(n+1)**3} tets={6*n**3} mode={mode} {stats}",flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared",type=Path)
    parser.add_argument("cli",type=Path)
    parser.add_argument("--reference",type=Path)
    parser.add_argument("--time-exe",type=Path)
    parser.add_argument("--repeat",type=int,default=3)
    parser.add_argument("--sizes",type=int,nargs="+",default=[4,8,16])
    opts = parser.parse_args()
    if opts.repeat < 1 or any(n < 2 or n%2 for n in opts.sizes):
        parser.error("Use positive repetitions and even mesh subdivisions >= 2.")
    with tempfile.TemporaryDirectory() as work:
        benchmark(opts,Path(work))
