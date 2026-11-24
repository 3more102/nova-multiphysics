# Forward polarization optics reference

The header-only C++20 API nova/physics/polarization.hpp models coherent transverse
Jones fields at normal incidence. Each layer has two passive complex principal
indices, a thickness in metres and an in-plane principal axis angle in radians.
The first index corresponds to the axis (cos(theta), sin(theta)); the caller
provides the indices (including any wavelength or director dependence).

With exp(-i omega t), the principal propagation factors are
a=exp(i k0 n_fast d) and b=exp(i k0 n_slow d). Each layer is
J=R(theta) diag(a,b) R(theta)^T. The total is J_last ... J_first.
apply_jones evaluates this matrix on a complex transverse field. For equal
reference impedances, the squared field norm ratio gives relative transmitted
intensity in this forward-only model. Each passive layer is a contraction
in the Euclidean Jones norm. Lossless layers preserve that norm.

This is a reflection-free reference model, not an anisotropic multilayer
Maxwell scattering solver. It excludes impedance discontinuities, multiple
reflections, oblique propagation, out-of-plane director effects, scattering,
depolarization and nonlinear self-consistency. It should not be used as an
absolute device transmission prediction where these effects are significant.

Verification: empty cascade, isotropic phase, 45 degree half- and quarter-wave
plates, Beer-Lambert attenuation, noncommuting cascade order, gain and NaN
rejection. CTest target nova_polarization_reference_tests.
