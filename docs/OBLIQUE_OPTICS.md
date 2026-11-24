# Oblique-incidence optics (experimental reference)
 
The header `nova/physics/oblique_optics.hpp` extends the existing normal-incidence reference models with coherent oblique TE and TM scattering in planar isotropic nonmagnetic multilayers. Wavelength and thickness are in metres, angle in radians from the normal, and exterior refractive indices are real positive. Complex passive layers use exp(-i omega t).
 
For conserved normalized tangential wavenumber s=n_inc sin(theta), the outgoing longitudinal wavenumber is k0 q, with q=sqrt(n^2-s^2) on the passive branch. Normalized tangential electric-field admittance is Y=q for TE and Y=n^2/q for TM. Fresnel amplitudes use r=(Y_a-Y_b)/(Y_a+Y_b), t=2Y_a/(Y_a+Y_b). The recursive cascade uses only exp(i k0 q d), avoiding exponentially growing transfer entries in thick absorbing or evanescent layers. TM amplitudes are tangential electric-field ratios, not the convention using oppositely oriented reflected p unit vectors.
 
R=|r|^2, T=Re(Y_exit)/Re(Y_inc)|t|^2 for a propagating exit, T=0 for an evanescent exit, and A=1-R-T. A is a residual, not independently integrated absorption. Small negative rounding residuals are not clipped. Exact critical/grazing longitudinal wavevectors, singular intermediate interfaces and numerically invalid power balances are rejected. Near-critical or resonant conditioning requires independent checking. No anisotropy, gain, magnetic materials, diffraction, finite beams, local field recovery or full-wave FEM is implemented.
 
Verification target: nova_oblique_optics_tests. Cases cover normal-incidence agreement with the existing implementation, Brewster angle, TE Fresnel reflection, total internal reflection, analytical frustrated TIR, thick absorbing metal, zero-thickness invariance, reciprocal transmission, and invalid input.
 
Run after configuring and building:
 
```sh
ctest --test-dir build -R nova_oblique_optics_tests --output-on-failure
```
 
Potential applications include angular DBR/filter studies, prism-coupled planar metal films, and reference calculations for future Maxwell FEM. These applications are not a validated general plasmonic or waveguide solver.
