# Passive optical dispersion (experimental)

The installed C++20 header `nova/physics/optical_dispersion.hpp` supplies isotropic,
nonmagnetic frequency-dependent material values for the existing driven 3D Maxwell
solver. It is header-only; the existing install rule exports it.

The time convention is exp(-i omega t). All angular frequencies are in rad/s,
vacuum wavelengths are in metres, and relative permittivity is dimensionless.

    epsilon_r = epsilon_inf
                - sum(wp^2 / (omega^2 + i gamma omega))
                + sum(f w0^2 / (w0^2 - omega^2 - i gamma omega))

Here wp is plasma frequency, w0 resonance frequency, f nonnegative dimensionless
oscillator strength, and gamma nonnegative damping. The imaginary part is
nonnegative for passive media. Negative real permittivity is supported by the
driven Maxwell assembly, but not by the lossless positive-material cavity
eigenproblem.

Example (illustrative parameters, not a fit to any named metal):

```cpp
#include "nova/physics/optical_dispersion.hpp"
#include <array>

const double wavelength = 1.55e-6;
const std::array<nova::physics::DrudePole,1> poles{{{1.4e16,1.0e14}}};
const auto material = nova::physics::optical_drude_lorentz_material(
    wavelength, 1.0, poles);
// Pass this material per cell to assemble_maxwell_3d, with
// vacuum_wavenumber = 2*pi/wavelength.
```

Alternatively, optical_material_from_nk computes epsilon=(n+i k)^2.
interpolate_optical_index linearly interpolates nonnegative n and k separately
on a strictly increasing wavelength table. Extrapolation is rejected.
Interpolation preserves passivity but does not enforce Kramers-Kronig consistency.

Invalid or nonfinite parameters and undamped singular Lorentz poles are rejected.
Extreme floating-point ratios outside representable range are not supported.

The regression suite includes closed-form Drude and Lorentz examples, superposition,
wavelength conversion, table interpolation, invalid-input rejection and a
manufactured driven Maxwell solution with positive integrated absorption.

Scope limits: scalar local linear permittivity; no material database, parameter
fitting, anisotropy, nonlinear susceptibility, waveguide modes, ports, PML, or
frequency-dependent eigenmode solve. The material API alone does not validate
a real device or establish commercial readiness.
