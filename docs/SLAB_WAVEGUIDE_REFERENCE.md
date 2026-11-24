# Symmetric dielectric slab waveguide reference

The header-only C++20 API nova/physics/slab_waveguide.hpp calculates bound TE and TM modes of a symmetric, nonmagnetic, isotropic, lossless dielectric slab. Inputs are vacuum wavelength, full core thickness, core and cladding refractive indices. This is an independent semi-analytical reference, not a general FEM eigenmode solver.

For half-thickness a and k0=2*pi/lambda, define V=k0*a*sqrt(n_core^2-n_clad^2), u=a*h and w=a*alpha. The core transverse wavenumber h and cladding decay alpha satisfy u^2+w^2=V^2. The modal condition is u=m*pi/2+atan(q*w/u), where q=1 for TE and q=n_core^2/n_clad^2 for TM. The order m is a nonnegative integer; V>m*pi/2 is required for a bound mode.

The solver parameterizes u=V*cos(theta), w=V*sin(theta) and brackets the monotonic dispersion function on [0,pi/2]. It reports effective index, propagation constant beta, transverse wavenumber, evanescent decay constant, and absolute dimensionless equation residual. The effective index is hypot(n_clad,alpha/k0). Cutoff modes are excluded; numerically unresolvable modes or extreme normalized frequencies V>10000 throw rather than silently producing invalid results.

Verification: the test executable checks the exact TE0 u=w=pi/4 case, independent tangent/cotangent interface equations for both polarizations, dispersion identity, mode count and ordering, cutoff crossing, confinement trend and invalid inputs.

No claim is made for lossy, asymmetric, anisotropic, nonlinear, finite-cladding, leaky or arbitrary 2D cross-section guides. This benchmark is intended for independent verification of future vector eigenmode FEM implementations.
