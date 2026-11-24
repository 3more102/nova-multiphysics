# Passive optical dispersion

The header `nova/physics/optical_dispersion.hpp` evaluates a scalar
Drude-Lorentz relative permittivity for the exp(-i omega t) convention:

epsilon_r(omega) = epsilon_infinity +
sum_j wp_j^2 / (w0_j^2 - omega^2 - i gamma_j omega).

Each oscillator uses a nonnegative plasma frequency wp (square root of
oscillator strength), resonance frequency w0, and damping rate gamma,
all in rad/s. A zero resonance represents Drude free carriers.
Positive damping produces nonnegative imaginary permittivity at positive
frequency. The real part can be negative. Parameters must be supplied by
the caller; no measured gold or other material dataset is included.

Example at a vacuum wavelength of 1550 nm, with illustrative rather than
measured material parameters:

~~~cpp
#include "nova/physics/optical_dispersion.hpp"
#include "nova/physics/maxwell_3d.hpp"
#include <array>

using namespace nova::physics;
const double w = optical_angular_frequency(1550e-9);
const std::array<OpticalOscillator,1> metal{{{3*w,0.0,0.1*w}}};
const MaxwellCellMaterial3D m{optical_permittivity(w,1.0,metal),1.0};
// Supply m to the driven Maxwell assembler with k0 = w / 299792458.
~~~

Evaluate materials separately at each frequency in a driven frequency
sweep. The bounded lossless cavity eigenproblem does not accept complex,
negative or frequency-dependent eigenproblem coefficients; do not insert
a single evaluated dispersive permittivity into that solver and interpret
the result as a dispersive eigenfrequency.

Analytical tests cover negative Drude real permittivity, positive loss,
high-frequency response, resonant Lorentz loss, oscillator superposition,
scale invariance at large frequencies, rejection of exact undamped poles
and invalid parameters, SI wavelength conversion and driven Maxwell assembly.
Run the normal CMake/CTest suite. These checks are not a validation of
measured optical constants, a nonlinear optical model, a dispersive
eigensolver, or a complete plasmonic device workflow. No PML, ports,
anisotropy or gain model is introduced here.

The evaluation scales frequency products before multiplication and rejects
detected nonrepresentable terms. Extreme ratios may underflow and are not
guaranteed to be supported.
