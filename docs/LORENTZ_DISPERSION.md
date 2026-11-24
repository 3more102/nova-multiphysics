# Passive Lorentz dispersion reference

The header-only C++20 API `nova/physics/lorentz_dispersion.hpp` evaluates
isotropic relative permittivity and its analytic angular-frequency derivative.
For the exp(-i omega t) convention,

epsilon(omega) = epsilon_inf + sum S_j omega_j^2 /
(omega_j^2 - omega^2 - i gamma_j omega).

omega, omega_j and gamma_j are in rad/s; S_j is dimensionless.
epsilon_inf > 0, omega_j > 0, gamma_j >= 0 and S_j >= 0 are required.
For positive frequencies the imaginary part of each pole is nonnegative.
A nonzero undamped pole at its resonance is singular and rejected.
Zero-strength poles contribute nothing.

The derivative is the sum of S_j omega_j^2 (2 omega + i gamma_j)
divided by the square of each complex denominator. The implementation
normalizes frequencies before squaring and rejects nonfinite outputs.
Extreme ratios may underflow in floating point; no arbitrary precision is claimed.

Use the returned relative_permittivity as the scalar complex material at
a single frequency in MaxwellCellMaterial3D (see MAXWELL_HCURL_3D.md).
The current real lossless cavity eigenproblem is not a dispersive eigenproblem.
This API does not fit experimental data, solve auxiliary-field time-domain
equations or implement nonlinear self-consistency.

Tests cover exact static and resonant responses and derivatives, high-frequency
limit, multiple poles at optical-scale frequencies, finite-difference
derivative comparison, passivity and invalid/singular input.
Run the existing NOVA CMake/CTest suite. CI results must be checked before
merging; these are mathematical reference checks, not experimental validation.
