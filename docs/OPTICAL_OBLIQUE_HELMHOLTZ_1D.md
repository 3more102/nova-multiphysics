# Oblique-incidence optical Helmholtz FEM

The C++20 solve_oblique_helmholtz_1d API extends the normal-incidence scalar solver to planar isotropic, nonmagnetic media. Coordinates use s=k0*z, conserved transverse wavevector b=n_i*sin(theta), and exp(-i omega t). The outgoing longitudinal branch q=sqrt(epsilon-b^2) has nonnegative imaginary part and nonnegative real part for passive media.

The scalar unknown u is E_y (TE) or Z0*H_y (TM), independently normalized to unit incident scalar amplitude. It obeys -(a*u')'-c*u=0 with a=1 for TE and a=1/epsilon for TM, c=a*(epsilon-b^2). The consistent linear element matrix is a/h*[1,-1;-1,1]-c*h/6*[2,1;1,2]. Continuity of u and a*u' follows from the weak form. Exact outgoing boundary terms subtract i*g at each end, with g=q for TE and q/epsilon for TM; the left load is -2*i*g_i. The existing adjacent-pivoted complex tridiagonal solver preserves O(N) storage and arithmetic.

For scalar amplitudes r_s=u_left-1 and t_s=u_right, R=|r_s|^2 and T=Re(g_t)/g_i*|t_s|^2. The returned tangential electric amplitudes are r=r_s, t=t_s for TE and r=-r_s, t=(g_t/g_i)*t_s for TM, consistent with the analytical oblique stack API. Evanescent exit waves carry zero normal transmitted power.

Cell absorptance is computed independently from exact P1 integrals: TE uses Im(epsilon)*integral(|u|^2)/g_i; TM uses -Im(1/epsilon)*(integral(|u'|^2)+b^2*integral(|u|^2))/g_i. The sum A is not inferred from 1-R-T; the signed residual is reported separately. A small residual alone does not certify spatial accuracy.

Validation excludes gain, nonfinite inputs, grazing incidence, exact exit critical angle, zero TM permittivity, and invalid meshes. The solver is not a general vector Maxwell or PML solver; it does not support anisotropy or nonlinear self-consistency.

## Verification

The nova_optical_oblique_helmholtz_1d_tests CTest target compares both polarizations with analytical scattering for a lossy asymmetric three-layer stack, total internal reflection, Brewster incidence and normal incidence. It checks complex amplitudes, power fractions, passive cell loss, second-order refinement and invalid input.

An independent Python/SciPy banded FEM calculation compared with homogeneous-layer exact state propagation yielded the following transmission amplitude errors for the lossy stack:

| Elements per layer | TE | TM |
|---|---:|---:|
| 200 | 3.87028e-5 | 3.18781e-5 |
| 400 | 9.67596e-6 | 7.96974e-6 |

Across the eight polarization/geometry cases at 400 elements per layer, the maximum amplitude or power fraction discrepancy was 9.676e-6 and the maximum absolute power balance residual was 3.360e-12. These are independent numerical checks, not a claim that the repository C++ CI has passed.

## Frustrated total internal reflection

A finite lossless low-index barrier between identical high-index propagating ports supports evanescent tunneling. Set exterior index n=1.5, barrier index m=1, angle theta=0.95 rad, b=n*sin(theta), q=sqrt(n^2-b^2), and kappa=sqrt(b^2-m^2). Define g=q and v=kappa for TE, or g=q/n^2 and v=kappa/m^2 for TM. Exact interface matching gives

    T = 1 / (1 + ((g^2+v^2)/(2*g*v))^2 * sinh(2*pi*kappa*d/lambda)^2)
    R = 1 - T, A = 0.

This reference does not call the analytical multilayer implementation. Regression cases cover both polarizations at four normalized thicknesses, positive and decreasing transmitted power, lossless energy conservation, and mesh convergence between 100 and 400 elements.

Independent SciPy banded P1 calculation (not the repository C++ executable):

| d/lambda | Exact TE T | Exact TM T |
|---|---:|---:|
| 0.05 | 0.9510579034 | 0.9362439687 |
| 0.25 | 0.3491434231 | 0.2884484956 |
| 0.75 | 0.0052291142 | 0.0039566412 |
| 1.5 | 0.00000721304 | 0.00000545083 |

The largest absolute transmission error across the eight cases at 400 elements was 2.453e-7; all independently computed power balance residuals were below 1e-10. This benchmark distinguishes lossless evanescent decay from actual absorption and from an evanescent exterior with zero normal transmitted power.
