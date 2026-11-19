"""Verify tensor heat assembly, step-doubling decisions and final temperatures.

Run: python verification/compare_adaptive_heat_3d.py build/nova_anisotropic_heat_3d
"""
import argparse
import csv
import io
import subprocess

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import factorized, spsolve
from skfem import Basis, BilinearForm, ElementTetP1, LinearForm, asm
from skfem.models.poisson import mass

from compare_transient_heat_3d import kuhn_cube


END_TIME = 0.02
TENSOR = np.array([[2.0, 0.4, -0.25], [0.4, 1.5, 0.3], [-0.25, 0.3, 1.0]])


@BilinearForm
def conduction(u, v, w):
    return (1 + 0.5*w.x[0]) * np.einsum("i...,ij,j...->...", v.grad, TENSOR, u.grad)


@LinearForm
def volume_source(v, w):
    return (1 + w.x[2]) * v


def independent_system(n):
    coords, mesh = kuhn_cube(n)
    basis = Basis(mesh, ElementTetP1(), intorder=2)
    free = basis.complement_dofs(basis.get_dofs())
    m = asm(mass, basis)[free][:, free].tocsc()
    k = asm(conduction, basis)[free][:, free].tocsc()
    f = asm(volume_source, basis)[free]
    xyz = coords[free]
    initial = 64 * np.prod(xyz * (1 - xyz), axis=1) * (1 + 0.5*xyz[:, 0])
    return coords, free, m, k, f, initial


def radau_reference(m, k, f, initial):
    # The independent integrator chooses its own time grid. The dense
    # Jacobian is limited to 27 free DOFs in this verification fixture.
    mass_solve = factorized(m)
    result = solve_ivp(lambda t, u: mass_solve((1+t)*f - k@u), (0, END_TIME),
                       initial, method="Radau", jac=-mass_solve(k.toarray()),
                       rtol=1e-11, atol=1e-13, t_eval=[END_TIME])
    if not result.success:
        raise RuntimeError(result.message)
    return result.y[:, -1]


def read_trial_run(executable, scheme, n, tolerance):
    output = subprocess.check_output([executable, scheme, str(n), str(tolerance), "--csv"],
                                     text=True)
    trials = []
    rows = csv.reader(io.StringIO(output))
    for row in rows:
        if row == ["x", "y", "z", "temperature"]:
            break
        if len(row) != 5 or row[0] != "#trial" or row[4] not in ("0", "1"):
            raise ValueError(f"Unexpected trial record {row}")
        trials.append((float(row[1]), float(row[2]), float(row[3]), row[4] == "1"))
    field = np.array(list(rows), dtype=float)
    return trials, field


def check_run(executable, n, scheme, tolerance, system, reference):
    coords, free, m, k, f, initial = system
    trials, field = read_trial_run(executable, scheme, n, tolerance)
    np.testing.assert_allclose(field[:, :3], coords, rtol=0, atol=2e-15)
    theta = 1.0 if scheme == "be" else 0.5
    divisor = 1 if scheme == "be" else 3

    def step(state, t, dt):
        rhs = (m - (1-theta)*dt*k) @ state + dt*(1+t+theta*dt)*f
        return spsolve(m + theta*dt*k, rhs)

    time = 0.0
    state = initial.copy()
    accepted = 0
    for start, dt, nova_error, take_step in trials:
        np.testing.assert_allclose(start, time, rtol=0, atol=2e-15)
        if not (0 < dt and start + dt <= END_TIME + 2e-15):
            raise AssertionError("Trial outside the physical time interval")
        full = step(state, start, dt)
        half = step(state, start, dt/2)
        refined = step(half, start + dt/2, dt/2)
        # NOVA's norm includes the prescribed (zero-error) vertices.
        error = np.linalg.norm(refined-full) / (np.sqrt(len(coords)) * tolerance * divisor)
        np.testing.assert_allclose(nova_error, error, rtol=2e-7, atol=2e-8)
        if take_step != (error <= 1):
            raise AssertionError("Independent trial decision disagrees")
        if take_step:
            state = refined
            time = start + dt
            accepted += 1
    if accepted < 2 or accepted == len(trials):
        raise AssertionError("Fixture must exercise both rejection and acceptance")
    np.testing.assert_allclose(time, END_TIME, rtol=0, atol=2e-15)
    replay = np.zeros(len(coords))
    replay[free] = state
    np.testing.assert_allclose(field[:, 3], replay, rtol=0, atol=2e-11)
    discrepancy = np.max(np.abs(replay - field[:, 3]))
    temporal_error = np.max(np.abs(field[free, 3] - reference))
    print(f"n={n} {scheme} atol={tolerance:g} accepted={accepted} rejected={len(trials)-accepted} "
          f"max_replay_difference={discrepancy:.3e} max_Radau_error={temporal_error:.3e}")
    return temporal_error


def compare(executable):
    for n in (3, 4):
        system = independent_system(n)
        reference = radau_reference(*system[2:])
        for scheme, limit in (("be", 1e-3), ("cn", 5e-5)):
            loose, tight = [check_run(executable, n, scheme, tol, system, reference)
                            for tol in (1e-4, 1e-6)]
            if not (tight < loose/4 and tight < limit):
                raise AssertionError(f"{scheme} terminal error exceeds its fixture acceptance band")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable")
    compare(parser.parse_args().executable)
