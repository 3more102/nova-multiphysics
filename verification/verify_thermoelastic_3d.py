"""Check coupled stationary/transient Gmsh cases using exported VTK fields."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import math
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET


def run(exe, mesh, *args, error=None):
    proc = subprocess.run([str(exe), str(mesh), *map(str, args)], capture_output=True, text=True)
    if error:
        if proc.returncode == 0 or error not in proc.stderr:
            raise ValueError(f"Expected {error!r}, got {proc.returncode}: {proc.stderr}")
    elif proc.returncode:
        raise ValueError(f"{proc.args}\n{proc.stderr}")
    return proc.stdout


def read_vtk(path):
    # Count tokens by association and tuple width, including every tensor entry.
    lines = path.read_text(encoding="ascii").splitlines()
    if lines[0] != "# vtk DataFile Version 3.0" or lines[2:4] != ["ASCII", "DATASET UNSTRUCTURED_GRID"]:
        raise ValueError("Invalid VTK header")
    stream = iter(" ".join(lines[4:]).split())

    def expect(*words):
        for word in words:
            if next(stream) != str(word):
                raise ValueError(f"Malformed VTK record, expected {word!r}")

    def values(n):
        data = tuple(float(next(stream)) for _ in range(n))
        if not all(map(math.isfinite, data)):
            raise ValueError("Non-finite VTK field")
        return data

    expect("FIELD", "FieldData", 3)
    meta = {}
    for _ in range(3):
        name = next(stream)
        expect(1, 1, "double")
        meta[name] = values(1)[0]
    expect("POINTS")
    n = int(next(stream))
    expect("double")
    points = [values(3) for _ in range(n)]
    expect("CELLS")
    count = int(next(stream))
    expect(5 * count)
    cells = []
    for _ in range(count):
        expect(4)
        cell = tuple(int(next(stream)) for _ in range(4))
        if min(cell) < 0 or max(cell) >= n or len(set(cell)) != 4:
            raise ValueError("Invalid VTK connectivity")
        cells.append(cell)
    expect("CELL_TYPES", count)
    for _ in cells:
        expect(10)
    expect("POINT_DATA", n, "SCALARS", "temperature", "double", 1, "LOOKUP_TABLE", "default")
    temp = values(n)
    expect("VECTORS", "displacement", "double")
    disp = [values(3) for _ in points]
    expect("CELL_DATA", count, "SCALARS", "material_id", "int", 1, "LOOKUP_TABLE", "default")
    fields = {"tags": [int(next(stream)) for _ in cells]}
    token = next(stream)
    if token == "SCALARS":
        expect("heat_generation", "double", 1, "LOOKUP_TABLE", "default")
        fields["generation"] = values(count)
        token = next(stream)
    if token != "VECTORS":
        raise ValueError("Missing heat flux vectors")
    expect("heat_flux", "double")
    fields["flux"] = [values(3) for _ in cells]
    for kind, name, width in [("VECTORS", "temperature_gradient", 3),
                               ("TENSORS", "strain", 9), ("TENSORS", "mean_stress", 9)]:
        expect(kind, name, "double")
        fields[name] = [values(width) for _ in cells]
    if next(stream, None) is not None:
        raise ValueError("Unexpected trailing VTK data")
    return points, cells, temp, disp, fields, meta


def close(actual, expected, tolerance):
    if not math.isclose(actual, expected, abs_tol=tolerance, rel_tol=0):
        raise ValueError(f"Measured {actual:.16g}, expected {expected:.16g} within {tolerance:g}")


def diagonal(x, y, z):
    return (x, 0, 0, 0, y, 0, 0, 0, z)


def check_tetra(path, amplitude, conductivity, tag=42):
    points, cells, temp, disp, fields, meta = read_vtk(path)
    if len(cells) != 1 or fields["tags"] != [tag]:
        raise ValueError("Wrong Gmsh region mapping")
    s = 0.01 * amplitude
    for (_, _, z), t, u in zip(points, temp, disp):
        close(t, 300 + amplitude*z, 3e-10)
        for value, expected in zip(u, (0, 0, s*z/4)):
            close(value, expected, 1e-11)
    expected_fields = {
        "temperature_gradient": (0,0,amplitude),
        "flux": tuple(-k * amplitude for k in conductivity),
        "strain": diagonal(0,0,s/4),
        "mean_stress": diagonal(-100*s/4,-100*s/4,0),
    }
    for name, ref in expected_fields.items():
        for value, expected in zip(fields[name][0], ref):
            close(value, expected, 2e-9)
    # Integrate 0.5*E*[2*(s*z)^2+(s/4-s*z)^2] over 0<=x+y+z<=1.
    energy = 50*s*s*(3/60 - 0.5/24 + 1/96)
    close(meta["elastic_strain_energy_J"], energy, 3e-12)
    if meta["elastic_residual"] > 1e-9 or meta["elastic_iterations"] < 0:
        raise ValueError("Invalid mechanical solver diagnostics")
    return max(temp), meta["elastic_strain_energy_J"]


def verify(steady, transient, root, output):
    wall = root / "examples/gmsh_diffusion_3d/two_layer_wall.msh"
    tet = root / "examples/thermoelastic_3d/tetrahedron.msh"
    tet41 = root / "tests/data/unit_tetrahedron41.msh"
    law = ["--elastic","solid",100,0,0.01,300]
    fixture = ["--clamp","base",*law]
    # 1 + 2. Resistances in series with fixed cold temperature or film resistance.
    materials = ["--material",11,2,"--material",22,6]
    laws = ["--elastic",11,100,0,0.001,300,"--elastic",22,200,0.25,0.002,310]
    for mixed in (False, True):
        case = "mixed_wall" if mixed else "steady_wall"
        q = 240 if mixed else 300
        boundary = ["--convection",102,12,300] if mixed else ["--temperature",102,300]
        path = output / (case + ".vtk")
        run(steady,wall,*materials,"--temperature",101,400,*boundary,"--flux",0,0,
            *laws,"--clamp",0,"--vtk",path,"--balance-report")
        points,cells,temp,disp,fields,meta = read_vtk(path)
        for (x,_,_), t, u in zip(points,temp,disp):
            close(t, 400-q*(min(x,0.5)/2+max(0,x-0.5)/6), 3e-10)
            for v in u:
                close(v,0,1e-14)
        energy = 0.0
        max_stress_error = 0.0
        for e, cell in enumerate(cells):
            left = max(points[i][0] for i in cell) <= 0.5
            if fields["tags"][e] != (11 if left else 22):
                raise ValueError("Volume group mismatch")
            modulus,alpha,ref = (100,0.001,300) if left else (400,0.002,310)
            mean = sum(temp[i]-ref for i in cell)/4
            for v,w in zip(fields["mean_stress"][e],diagonal(*([-modulus*alpha*mean]*3))):
                close(v,w,1e-9)
                max_stress_error = max(max_stress_error,abs(v-w))
            for v,w in zip(fields["flux"][e],(q,0,0)):
                close(v,w,2e-8)
        # Cross section is one. Integrate (linear temperature - Tref)^2 per layer.
        tint = 400 - q/4
        tlast = tint - q/12
        for a,b,modulus,alpha,ref in [(400,tint,100,.001,300),(tint,tlast,400,.002,310)]:
            a,b = a-ref,b-ref
            energy += 1.5 * modulus * alpha**2 * .5 * (a*a+a*b+b*b)/3
        close(meta["elastic_strain_energy_J"],energy,3e-10)
        print(f"{case} Tinterface={tint:g} q={fields['flux'][0][0]:.12g} "
              f"energy={meta['elastic_strain_energy_J']:.12g} max_stress_error={max_stress_error:.3g}")

    # 3. Exact T=300+10z for an SPD tensor and nonzero Neumann data on three faces.
    aniso = ["--material-tensor","solid",2,3,4,.2,.5,-.25,"--temperature","base",300,
             "--flux","x_min",5,"--flux","y_min",-2.5,"--flux","incline",-42.5/math.sqrt(3)]
    path = output/"anisotropic.vtk"
    run(steady,tet,*aniso,*fixture,"--vtk",path)
    print("anisotropic Tpeak and elastic energy",check_tetra(path,10,(.5,-.25,4)))

    # 4. A cooling body with one free thermal DOF. The reference P1 element has
    # M44=C/60, K44=k/6, so the amplitude solves a'=-10a. Both schemes have
    # a known amplification factor. Check the physical cooling and every FEM frame.
    for scheme, ratio, mesh, group, tag in (
        ("be",1/(1+.25),tet,"base",42),
        ("cn",(1-.125)/(1+.125),tet41,"fixed-temperature",11)
    ):
        args = ["--material","solid",1,"--capacity","solid",1,"--initial",310,
                "--temperature",group,300,"--dt",.025,"--steps",8,"--scheme",scheme,
                "--clamp",group,*law,"--elastic-solver","ic0",1e-11,100]
        collection = output/f"cooling_{scheme}.pvd"
        budget = output/f"cooling_{scheme}.csv"
        run(transient,mesh,*args,"--vtk-series",collection,"--energy-csv",budget)
        datasets = ET.parse(collection).findall("./Collection/DataSet")
        if len(datasets) != 9:
            raise ValueError("Missing accepted thermal snapshot")
        rows = list(csv.DictReader(budget.open(encoding="ascii")))
        if len(rows) != 8:
            raise ValueError("Thermal energy row count mismatch")
        for step, ds in enumerate(datasets):
            close(float(ds.get("timestep")), .025*step,1e-14)
            amplitude = 10*ratio**step
            peak, energy = check_tetra(output/ds.attrib["file"],amplitude,(0,0,1),tag)
            if step:
                row = rows[step-1]
                close(float(row["heat_content_J"]), 300/6+amplitude/24, 3e-10)
                close(float(row["balance_residual_J"]),0,3e-10)
                close(float(row["dirichlet_energy_J"]),amplitude*(1-1/ratio)/24,3e-10)
        print(f"{scheme} t=0.2 Tpeak={peak:.12g} continuum_time_error={abs(peak-300-10*math.exp(-2)):.6g} "
              f"elastic_energy={energy:.12g} max_heat_defect={max(abs(float(r['balance_residual_J'])) for r in rows):.4g}")
        final = output/f"cooling_{scheme}_final.vtk"
        run(transient,mesh,*args,"--vtk",final)
        check_tetra(final,10*ratio**8,(0,0,1),tag)

        errors = [abs(peak-300-10*math.exp(-2))]
        with tempfile.TemporaryDirectory(dir=output) as refinement:
            for steps in (16,32):
                refined = args.copy()
                refined[refined.index("--dt")+1] = .2/steps
                refined[refined.index("--steps")+1] = steps
                path = Path(refinement)/f"refine{steps}.vtk"
                run(transient,mesh,*refined,"--vtk",path)
                points,_,temps,disp,_,_ = read_vtk(path)
                tip = next(i for i,p in enumerate(points) if p[2] == 1)
                errors.append(abs(temps[tip]-300-10*math.exp(-2)))
                close(disp[tip][2],.0025*(temps[tip]-300),1e-12)
        orders = [math.log2(a/b) for a,b in zip(errors,errors[1:])]
        if min(orders) < (0.9 if scheme == "be" else 1.95):
            raise ValueError(f"Insufficient {scheme} temporal order: {orders}")
        print(f"{scheme} time_refinement_errors={errors} orders={orders}")

    # Three symmetry planes leave free thermal expansion and fix rigid motion.
    free = ["--fix-component","x_min","x",0,"--fix-component","y_min","y",0,
            "--fix-component","base","z",0]
    path = output/"free_expansion.vtk"
    run(steady,tet,"--material","solid",1,"--temperature","base",305,*law,*free,"--vtk",path)
    points,_,temps,disps,fields,meta = read_vtk(path)
    for point,t,u in zip(points,temps,disps):
        close(t,305,1e-11)
        for coord,value in zip(point,u):
            close(value,.05*coord,1e-12)
    for value in fields["mean_stress"][0]:
        close(value,0,1e-10)
    close(meta["elastic_strain_energy_J"],0,1e-20)
    print(f"free expansion max_stress={max(map(abs,fields['mean_stress'][0])):.3g} Pa")

    prior = (output/"anisotropic.vtk").read_bytes()
    cases = [(law,"rigid modes"),(["--clamp","base"],"Missing elastic"),
             ([*fixture,*law],"Duplicate elastic"),([*fixture,"--clamp",777],"Unknown clamp"),
             (["--clamp","base","--elastic","solid",100,.5,.01,300],"isotropic elasticity"),
             ([*fixture,"--elastic-solver","jacobi",.1,0],"iteration limit")]
    cases.extend([([*law,*free[:-4]],"rigid modes"),
                  ([*law,*free,"--fix-component","base","w",0],"must be x, y or z"),
                  ([*fixture,"--fix-component","x_min","x",1],"Conflicting prescribed")])
    for extra, error in cases:
        run(steady,tet,*aniso,*extra,"--vtk",output/"anisotropic.vtk",error=error)
    run(steady,tet,"--material","solid",1,*fixture,"--vtk",output/"anisotropic.vtk",error="Absolute temperature")
    run(steady,wall,*materials,"--temperature",101,400,"--clamp",101,*laws,
        "--elastic-solver","jacobi",1e-12,1,"--vtk",output/"anisotropic.vtk",error="Thermoelastic CG failed")
    if prior != (output/"anisotropic.vtk").read_bytes():
        raise ValueError("Failed model changed a previously exported solution")

    # Competing runs with different temperatures must publish exactly one
    # coherent series. Check the amplitude of every frame against the winner.
    with tempfile.TemporaryDirectory(dir=output) as work:
        series = Path(work)/"race.pvd"
        def compete(temp):
            cmd = [str(transient),str(tet),"--material","solid","1","--capacity","solid","1",
                   "--initial",str(temp),"--temperature","base","300","--dt",".025","--steps","3",
                   *map(str,fixture),"--vtk-series",str(series)]
            return subprocess.run(cmd,capture_output=True).returncode
        starts = [301,306,312,320]
        with ThreadPoolExecutor(max_workers=len(starts)) as executor:
            winners = [t for t,rc in zip(starts,executor.map(compete,starts)) if rc == 0]
        if len(winners) != 1:
            raise ValueError(f"Expected one publishing run, found {winners}")
        frames = ET.parse(series).findall("./Collection/DataSet")
        if len(frames) != 4:
            raise ValueError("Race produced an incomplete PVD")
        for i,frame in enumerate(frames):
            check_tetra(Path(work)/frame.attrib["file"],(winners[0]-300)*.8**i,(0,0,1))



if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("steady",type=Path)
    cli.add_argument("transient",type=Path)
    cli.add_argument("--output",type=Path,help="Keep cases in a new or empty directory")
    args = cli.parse_args()
    repo = Path(__file__).resolve().parents[1]
    if args.output:
        args.output.mkdir(parents=True,exist_ok=True)
        if any(args.output.iterdir()):
            cli.error("--output must be empty")
        verify(args.steady.resolve(),args.transient.resolve(),repo,args.output)
    else:
        with tempfile.TemporaryDirectory() as work:
            verify(args.steady.resolve(),args.transient.resolve(),repo,Path(work))
