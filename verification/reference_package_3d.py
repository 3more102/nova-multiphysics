"""Independent double-precision DC reference for the supplied package at its final NOVA temperature."""
import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def compare(frame, terminal_current):
    points = np.asarray(frame["points"])
    cells = np.asarray(frame["cells"])
    tags = np.asarray(frame["material_id"], dtype=int).ravel()
    temperature = np.asarray(frame["temperature"]).ravel()
    v_nova = np.asarray(frame["electric_potential"]).ravel()
    # Independent constitutive evaluation for this specific reproducible coupon.
    sigma = np.zeros(len(cells))
    mean = temperature[cells].mean(axis=1)
    for tag, ref, slope in [(1, 5.8e7, -.0039), (4, 1000., -.003), (5, 5.8e7, -.0039)]:
        idx = tags == tag
        sigma[idx] = ref*np.exp(slope*(mean[idx]-300))
    active = np.unique(cells[sigma > 0])
    xyz = points[cells]
    affine = np.concatenate([np.ones((*xyz.shape[:2], 1)), xyz], axis=2)
    # Columns of inv(affine) contain the four physical shape polynomials.
    grads = np.linalg.inv(affine)[:, 1:, :].transpose(0, 2, 1)
    volume = np.abs(np.linalg.det(xyz[:, 1:]-xyz[:, :1]))/6
    local = (sigma*volume)[:, None, None] * (grads @ grads.transpose(0, 2, 1))
    rows = np.repeat(cells, 4, axis=1).ravel()
    cols = np.tile(cells, (1, 4)).ravel()
    A = coo_matrix((local.ravel(), (rows, cols)), shape=(len(points),)*2).tocsr()
    zmax = points[:, 2].max()
    high = (points[:, 2] == zmax) & (points[:, 0] >= .003-1e-14)
    low = ((points[:, 2] == zmax) & (points[:, 0] <= .001+1e-14)) | (points[:, 2] == 0)
    free = np.setdiff1d(active, np.flatnonzero(high | low))
    V = np.zeros(len(points))
    V[high] = .75
    V[free] = spsolve(A[free][:, free], -(A @ V)[free])
    voltage_error = np.abs(V[active]-v_nova[active]).max()
    current = float((A @ V)[high].sum())
    gradient = np.einsum("ei,eij->ej", V[cells], grads)
    joule = np.einsum("ij,ij->i", gradient, gradient)*sigma
    power = float(volume @ joule)
    power_nova = float(volume @ np.asarray(frame["joule_density"]).ravel())
    assert voltage_error < 1e-8
    assert abs(current-terminal_current) < 2e-7 and abs(power-power_nova) < 1e-8
    assert math.isclose(current*.75, power, abs_tol=1e-8)
    return dict(voltage_max_error_V=float(voltage_error), terminal_current_A=current,
                terminal_current_error_A=abs(current-terminal_current), joule_W=power,
                joule_error_W=abs(power-power_nova))


def mechanical_and_transient(frame):
    """Independently assemble small-strain mechanics and a continuous-time heat/DC DAE."""
    from itertools import combinations
    from scipy.integrate import solve_ivp
    points = np.asarray(frame["points"])
    cells = np.asarray(frame["cells"])
    tags = np.asarray(frame["material_id"], dtype=int).ravel()-1
    xyz = points[cells]
    ne, n = len(cells), len(points)
    shape = np.concatenate([np.ones((ne, 4, 1)), xyz], axis=2)
    grad = np.linalg.inv(shape)[:, 1:, :].transpose(0, 2, 1)
    volume = np.abs(np.linalg.det(xyz[:, 1:]-xyz[:, :1]))/6
    temperature = np.asarray(frame["temperature"]).ravel()
    delta = temperature[cells].mean(axis=1)-300
    E = np.take([110e9, 2e9, 300e9, 130e9, 110e9, 4e9], tags)
    E *= np.exp(np.take([-.0004, -.001, -.0002, -.0003, -.0004, -.001], tags)*delta)
    nu = np.take([.34, .35, .22, .28, .34, .30], tags)
    alpha = np.take([16.5e-6, 30e-6, 7e-6, 2.6e-6, 16.5e-6, 20e-6], tags)
    lame = E*nu/((1+nu)*(1-2*nu))
    mu = E/(2*(1+nu))
    D = np.zeros((ne, 6, 6))
    D[:, :3, :3] = lame[:, None, None]
    for i in range(6):
        D[:, i, i] += (2 if i < 3 else 1)*mu
    # Engineering shear ordering xx yy zz xy yz xz.
    B = np.zeros((ne, 6, 12))
    for i in range(3):
        B[:, i, i::3] = grad[:, :, i]
    for a, (i, j) in enumerate([(0, 1), (1, 2), (0, 2)], 3):
        B[:, a, i::3], B[:, a, j::3] = grad[:, :, j], grad[:, :, i]
    dofs = (3*cells[:, :, None] + np.arange(3)).reshape(ne, 12)
    K = volume[:, None, None]*(B.transpose(0, 2, 1) @ D @ B)
    A = coo_matrix((K.ravel(), (np.repeat(dofs, 12, axis=1).ravel(), np.tile(dofs, (1, 12)).ravel())),
                   shape=(3*n, 3*n)).tocsr()
    thermal = np.zeros((ne, 6))
    thermal[:, :3] = (alpha*delta)[:, None]
    loads = np.einsum("e,eij,ejk,ek->ei", volume, B.transpose(0, 2, 1), D, thermal)
    rhs = np.bincount(dofs.ravel(), weights=loads.ravel(), minlength=3*n)
    free_u = np.repeat(points[:, 2] != 0, 3)
    u = np.zeros(3*n)
    u[free_u] = spsolve(A[free_u][:, free_u], rhs[free_u])
    stress = np.einsum("eij,ej->ei", D, np.einsum("eij,ej->ei", B, u[dofs])-thermal)
    disp_error = float(np.max(np.abs(u.reshape(n, 3)-frame["displacement"])))
    stress_error = float(np.max(np.abs(stress-np.asarray(frame["mean_stress"])[:, [0, 4, 8, 1, 5, 2]])))
    reaction_error = float(np.max(np.abs((A@u-rhs).reshape(n, 3)-frame["support_reaction"])))
    assert disp_error < 1e-15 and stress_error < .01 and reaction_error < 1e-8

    # Rebuild all thermal and electrical coefficients from the reference case.
    ref_sig = np.take([5.8e7, 0, 0, 1000., 5.8e7, 0], tags)
    b_sig = np.take([-.0039, 0, 0, -.003, -.0039, 0], tags)
    k0 = np.take([390., 3., 24., 130., 390., .8], tags)
    bk = np.take([-.0005, 0, -.001, -.002, -.0005, 0], tags)
    C0 = np.take([3.45e6, 1.8e6, 3e6, 1.65e6, 3.45e6, 1.5e6], tags)
    bc = np.take([.0003, .0008, .0004, .0005, .0003, .001], tags)
    G = volume[:, None, None]*(grad @ grad.transpose(0, 2, 1))
    M = volume[:, None, None]*(np.ones((4, 4))+np.eye(4))/20
    rows, cols = np.repeat(cells, 4, axis=1).ravel(), np.tile(cells, (1, 4)).ravel()

    def sparse(local):
        return coo_matrix((local.ravel(), (rows, cols)), shape=(n, n)).tocsr()

    zmax = points[:, 2].max()
    high = (points[:, 2] == zmax) & (points[:, 0] >= .003-1e-14)
    low = (points[:, 2] == 0) | ((points[:, 2] == zmax) & (points[:, 0] <= .001+1e-14))
    active = np.unique(cells[ref_sig > 0])
    free = np.setdiff1d(active, np.flatnonzero(high | low))
    faces = np.vstack([cells[:, list(c)] for c in combinations(range(4), 3)])
    unique, counts = np.unique(np.sort(faces, axis=1), axis=0, return_counts=True)
    surf = unique[counts == 1]
    base = np.all(points[surf, 2] == 0, axis=1)
    lid = (np.all(points[surf, 2] == zmax, axis=1)
           & (points[surf, 0].mean(axis=1) > .001) & (points[surf, 0].mean(axis=1) < .003))
    faces = surf[base | lid]
    htc = np.where(base[base | lid], 2500., 12.)
    eps = np.where(base[base | lid], 0., .7)
    a, b, c = points[faces].transpose(1, 0, 2)
    area = .5*np.linalg.norm(np.cross(b-a, c-a), axis=1)
    gp, gw = np.polynomial.legendre.leggauss(4)
    gp, gw = (1+gp)/2, gw/2
    N = np.array([[1-x, x*(1-y), x*y] for x in gp for y in gp])
    W = np.array([2*x*wx*wy for x, wx in zip(gp, gw) for wy in gw])

    def ode(t, excess):
        d = excess[cells]
        mean = d.mean(axis=1)
        sigma = ref_sig*np.exp(b_sig*mean)
        K = sparse(sigma[:, None, None]*G)
        v = np.zeros(n)
        v[high] = np.interp(t, [0, .01, .05, .06], [0, 1.5, 1.5, .75])
        v[free] = spsolve(K[free][:, free], -(K@v)[free])
        dv = np.einsum("en,enj->ej", v[cells], grad)
        q = volume*sigma*np.einsum("ei,ei->e", dv, dv)
        cond = (k0*np.exp(bk*mean))[:, None]*np.einsum("eij,ej->ei", G, d)
        rhs = np.bincount(cells.ravel(), weights=(q[:, None]/4-cond).ravel(), minlength=n)
        diff = excess[faces] @ N.T
        flux = htc[:, None]*diff + eps[:, None]*5.670374419e-8*((300+diff)**4-300.**4)
        cooling = area[:, None]*((flux*W) @ N)
        np.add.at(rhs, faces, -cooling)
        mass = sparse(M*(C0[:, None]*np.exp(bc[:, None]*d))[:, None, :])
        return spsolve(mass, rhs)

    def integrate(rtol, atol):
        x, calls = np.zeros(n), 0
        for t0, t1 in zip([0, .01, .05, .06], [.01, .05, .06, .08]):
            result = solve_ivp(ode, (t0, t1), x, method="BDF", atol=atol, rtol=rtol)
            assert result.success
            x, calls = result.y[:, -1], calls+result.nfev
        return 300+x, calls

    standard, _ = integrate(3e-6, 3e-8)
    refined, calls = integrate(3e-8, 3e-10)
    reference_change = float(np.max(np.abs(standard-refined)))
    thermal_error = float(np.max(np.abs(temperature-refined)))
    assert reference_change < 1e-4 and thermal_error < .05
    return dict(displacement_max_error_m=disp_error, stress_max_error_Pa=stress_error,
                reaction_max_error_N=reaction_error, final_temperature_max_error_K=thermal_error,
                reference_tightening_change_K=reference_change, reference_Tmax_K=float(max(refined)),
                refined_ode_nfev=calls)
