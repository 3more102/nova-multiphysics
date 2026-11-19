"""Compare every cooling-cube temperature with independent scikit-fem assembly and SciPy LU.

Run: python verification/compare_transient_heat_3d.py build/nova_transient_heat_3d
Reference environment: numpy, scipy and scikit-fem (optional verification tools).
"""
import argparse
import io
import itertools
import subprocess

import numpy as np
import scipy
from scipy.sparse.linalg import factorized
import skfem
from skfem import Basis, ElementTetP1, MeshTet, asm
from skfem.models.poisson import laplace, mass


def kuhn_cube(n):
    # Independently construct Kuhn simplices as all monotone corner-to-corner paths.
    side = n + 1
    vertices = [(x / n, y / n, z / n) for z, y, x in itertools.product(range(side), repeat=3)]
    cells = []
    for z, y, x in itertools.product(range(n), repeat=3):
        origin = x + side * y + side**2 * z
        for a, b, c in itertools.permutations((1, side, side**2)):
            cells.append([origin, origin + a, origin + a + b, origin + a + b + c])
    coords = np.array(vertices)
    mesh = MeshTet(np.ascontiguousarray(coords.T), np.ascontiguousarray(np.array(cells).T))
    return coords, mesh


def compare(executable, n, steps):
    text = subprocess.check_output([executable, str(n), str(steps), "--csv"], text=True)
    field = np.loadtxt(io.StringIO(text), delimiter=",", skiprows=1)
    coords, mesh = kuhn_cube(n)
    np.testing.assert_allclose(field[:, :3], coords, rtol=0, atol=2e-15)
    basis = Basis(mesh, ElementTetP1(), intorder=2)
    free = basis.complement_dofs(basis.get_dofs())
    capacity = asm(mass, basis)[free][:, free]
    stiffness = asm(laplace, basis)[free][:, free]
    half_dt = 0.01 / steps
    step_solve = factorized((capacity + half_dt * stiffness).tocsc())
    right = capacity - half_dt * stiffness
    state = np.prod(np.sin(np.pi * coords[free]), axis=1)
    for _ in range(steps):
        state = step_solve(right @ state)
    reference = np.zeros(len(coords))
    reference[free] = state
    discrepancy = np.max(np.abs(field[:, 3] - reference))
    np.testing.assert_allclose(field[:, 3], reference, rtol=0, atol=1e-10)
    # Independent high-order quadrature measures error of the whole P1 field.
    quadrature = Basis(basis.mesh, ElementTetP1(), intorder=5)
    points = quadrature.global_coordinates()
    exact = np.exp(-3 * np.pi**2 * 0.02) * np.prod(np.sin(np.pi * points), axis=0)
    error = np.sqrt(np.sum((quadrature.interpolate(field[:, 3]) - exact)**2 * quadrature.dx))
    print(f"n={n:2d} dofs={len(coords):5d} tets={mesh.nelements:6d} "
          f"max_backend_difference={discrepancy:.4e} L2_error={error:.6e}")
    return error


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable")
    args = parser.parse_args()
    print(f"scikit-fem {skfem.__version__}, SciPy {scipy.__version__}, NumPy {np.__version__}")
    errors = [compare(args.executable, n, 20) for n in (3, 6, 12)]
    orders = np.log2(np.array(errors[:-1]) / errors[1:])
    print("L2 refinement orders:", *[f"{r:.4f}" for r in orders])
    if not np.all(orders > 1.6):
        raise SystemExit("Cooling-cube spatial refinement below the regression threshold.")
