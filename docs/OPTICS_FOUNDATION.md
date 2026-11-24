# Optics, photonics and plasmonics: initial reference module

The C++20 API `nova/physics/optics.hpp` provides analytical and semi-analytical reference models. It is not a Maxwell FEM or nonlinear wave solver.

## Conventions and scope
SI units, vacuum wavelength, nonmagnetic isotropic media, time dependence exp(-i omega t). Passive indices have nonnegative imaginary part. Exterior half spaces for multilayer propagation must be lossless and have positive real index. Layers may have complex passive refractive index, including purely imaginary index. Zero thickness is allowed.

## Linear optics and multilayer photonics
At normal incidence r_ab=(n_a-n_b)/(n_a+n_b), t_ab=1+r_ab. Cascade a downstream reflection r and transmission t using denominator 1+r_ab*r, then propagate through each finite layer with p=exp(i k0 n d). The recursive reflection accumulates p^2 and transmission p. This avoids exponentially growing transfer matrix entries for thick absorbing layers. R=|r|^2, T=(n_exit/n_incident)|t|^2, A=1-R-T (small negative rounding residuals may occur). Extremely ill-conditioned cases are not guaranteed. The normal-incidence API does not handle oblique incidence; see the extension below. Neither API handles anisotropy, diffraction, guided modes or field reconstruction.

## Material dispersion and plasmonics
Drude relative permittivity: epsilon_inf - omega_p^2/[omega(omega+i gamma)]. The passive square root gives the index. This is a parameterized model, not a fitted gold/silver database.
For real positive dielectric epsilon_d and passive metal Re(epsilon_m)<-epsilon_d, the planar TM SPP propagation constant is beta=k0 sqrt(epsilon_d epsilon_m/(epsilon_d+epsilon_m)). Decay constants are sqrt(beta^2-k0^2 epsilon_j) with positive real parts. Intensity propagation length is 1/(2 Im(beta)), infinite for a lossless mode. The model does not solve nanoparticles, gratings, arbitrary interfaces or leaky modes.

## Nonlinear optics
Prescribed-intensity scalar Kerr reference: delta_n=n2 I; delta_phi=k0 delta_n L. Both positive and negative n2 are supported. No self-consistent nonlinear field solution, frequency conversion, nonlinear absorption or pulse propagation.

## Verification
CTest target nova_optics_reference_tests checks Fresnel reflection, quarter-wave antireflection, the analytical quarter-wave DBR result, lossless power balance, thick absorbing layer stability, Drude signs, SPP dispersion and TM interface continuity, Kerr phase and invalid inputs. Run:
```sh
cmake --preset default
cmake --build --preset default
ctest --preset default --output-on-failure
```
Passing reference tests alone do not validate a general optical solver.

## Next implementation gates
1. Complex scalar Helmholtz FEM with manufactured solution and convergence tests, followed by absorbing boundary and port excitation.
2. Vector Maxwell H(curl) edge elements with cavity eigenfrequency benchmarks and spurious-mode checks.
3. Waveguide eigenmodes, effective index, confinement and loss; scattering ports, PML, dispersive material data and wavelength sweeps.
4. Self-consistent Kerr field iteration and chi(2)/chi(3) nonlinear models with continuation and convergence diagnostics.
5. Plasmonic waveguide and particle scattering, Mie and planar benchmarks, optical heating coupling; independent cross-solver comparisons.

## Oblique incidence extension
Include `nova/physics/oblique_optics.hpp` and call `oblique_incidence_stack` with a vacuum wavelength, positive real exterior indices, angle in radians, TE or TM polarization, and the existing Layer span. Reflection and transmission amplitudes refer to tangential electric field components (the reflected TM amplitude sign may differ from conventions using a local p-polarization basis).

For conserved transverse index s = n_incident sin(theta), use q = sqrt(n^2-s^2) on the outgoing/passive branch. Normalized admittances are Y_TE = q and Y_TM = n^2/q. The interface amplitude is (Y_a-Y_b)/(Y_a+Y_b); the existing scattering recursion uses exp(i k0 q d). Power transmission is Re(Y_exit)/Re(Y_incident) times |t|^2. Evanescent lossless exit media therefore have zero transmitted normal power; finite evanescent gaps can transmit via tunneling.

Limitations: planar isotropic nonmagnetic passive layers only; exact zero longitudinal wavenumber and grazing incidence are rejected. Ill-conditioned resonances are not guaranteed. This is not a full-wave FEM, PML, or port solver. The CTest target nova_oblique_optics_tests exercises both polarizations, Brewster reflection, total internal reflection, frustrated TIR, thick lossy media, quarter-wave matching, and invalid angle rejection.
