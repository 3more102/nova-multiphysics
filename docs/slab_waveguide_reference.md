# Symmetric dielectric slab reference

The C++20 API nova::physics::solve_symmetric_slab computes strictly guided
TE and TM planar modes in a symmetric lossless isotropic nonmagnetic slab.
It is an analytical dispersion reference for later FEM verification, not a
full-vector cross-sectional eigenmode solver, rib/ridge solver, or PML model.

Let a be the core half-thickness, k0 = 2 pi / wavelength, n1 > n2 > 0.
V = k0 a sqrt(n1^2 - n2^2); u^2 + w^2 = V^2.
Each mode m satisfies

    u = m pi/2 + atan2(f w, u)

where f = 1 for TE and f = n1^2/n2^2 for TM. Solve on the monotone
theta interval [0, pi/2], with u = V cos(theta), w = V sin(theta).
Only m pi/2 < V is guided. The cladding decay is w/a, and
n_eff^2 = n2^2 + (n1^2 - n2^2)(w/V)^2. Beta = k0 n_eff.

Even and odd scalar fields are cosine and sine respectively, continued
by exponentially decaying cladding fields. Their exact integrated core
and cladding longitudinal power fractions use the |E_y|^2 weight for TE
and |H_y|^2/epsilon for TM. Reported residuals are absolute phase residuals.
The reference intentionally bounds V to [1e-4, 1e4] and rejects
unrepresentable effective indices at cutoff rather than returning false modes.

Verification: independently computed SciPy brentq values at n1=1.5,
n2=1.0, V=2.5: TE0 1.4153986308765596, TE1 1.160423010582448,
TM0 1.3817754032255871, TM1 1.085752934406741.
C++ tests check these values, dispersion residual, mode cutoff,
scale invariance, physical confinement, and invalid inputs.
