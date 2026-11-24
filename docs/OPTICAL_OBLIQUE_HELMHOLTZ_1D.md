# Oblique-incidence optical Helmholtz FEM

The C++20 solve_oblique_helmholtz_1d API extends the normal-incidence scalar solver to planar isotropic, nonmagnetic media. Coordinates use s=k0*z, conserved transverse wavevector b=n_i*sin(theta), and exp(-i omega t). The outgoing longitudinal branch q=sqrt(epsilon-b^2) has nonnegative imaginary part and nonnegative real part for passive media.

The scalar unknown u is E_y (TE) or Z0*H_y (TM), independently normalized to unit incident scalar amplitude. It obeys -(a*u')'-c*u=0 with a=1 for TE and a=1/epsilon for TM, c=a*(epsilon-b^2). The consistent linear element matrix is a/h*[1,-1;-1,1]-c*h/6*[2,1;1,2]. Continuity of u and a*u' follows from the weak form. Exact outgoing boundary terms subtract i*g at each end, with g=q for TE and q/epsilon for TM; the left load is -2*i*g_i. The existing adjacent-pivoted complex tridiagonal solver preserves O(N) storage and arithmetic.

For scalar amplitudes r_s=u_left-1 and t_s=u_right, R=|r_s|^2 and T=Re(g_t)/g_i*|t_s|^2. The returned tangential electric amplitudes are r=r_s, t=t_s for TE and r=-r_s, t=(g_t/g_i)*t_s for TM, consistent with the analytical oblique stack API. Evanescent exit waves carry zero normal transmitted power.

Cell absorptance is computed independently from exact P1 integrals: TE uses Im(epsilon)*integral(|u|^2)/g_i; TM uses -Im(1/epsilon)*(integral(|u'|^2)+b^2*integral(|u|^2))/g_i. The sum A is not inferred from 1-R-T; the signed residual is reported separately. A small residual alone does not certify spatial accuracy.

Validation excludes gain, nonfinite inputs, grazing incidence, zero TM permittivity, and invalid meshes. The solver is not a general vector Maxwell or PML solver; it does not support anisotropy or nonlinear self-consistency.

At the exact exit critical angle, q_t=0 and the limiting outgoing boundary condition is homogeneous Neumann (zero scalar flux). Both polarizations give T=0. TE may have finite tangential E at the exit; TM has zero tangential E there, despite a potentially finite magnetic scalar field. The analytical stack API retains its separate critical-angle restriction.

The critical-exit regression uses n_i=1.5, theta=0.7, n_t=n_i*sin(theta), and a slab of thickness 0.37 wavelengths with 80 elements. With epsilon=n_t^2, the exact scalar field is u=2 everywhere, R=1 and A=0. A lossy slab (epsilon=2.1+0.3i) tests independent absorption and continuity from both sides of the critical angle (theta +/- 1e-9). An independent SciPy banded P1 computation gave A=0.6275856269 (TE) and A=0.7341033682 (TM), maximum exact-critical energy residual 4.17e-14 and maximum near-critical amplitude/absorption deviation 2.34e-4. These are numerical reference checks, not C++ CI results.

## Verification

The nova_optical_oblique_helmholtz_1d_tests CTest target compares both polarizations with analytical scattering for a lossy asymmetric three-layer stack, total internal reflection, Brewster incidence and normal incidence. It checks complex amplitudes, power fractions, passive cell loss, second-order refinement and invalid input.

An independent Python/SciPy banded FEM calculation compared with homogeneous-layer exact state propagation yielded the following transmission amplitude errors for the lossy stack:

| Elements per layer | TE | TM |
|---|---:|---:|
| 200 | 3.87028e-5 | 3.18781e-5 |
| 400 | 9.67596e-6 | 7.96974e-6 |

Across the eight polarization/geometry cases at 400 elements per layer, the maximum amplitude or power fraction discrepancy was 9.676e-6 and the maximum absolute power balance residual was 3.360e-12. These are independent numerical checks, not a claim that the repository C++ CI has passed.
