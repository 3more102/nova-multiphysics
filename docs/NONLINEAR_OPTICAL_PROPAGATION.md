# Scalar nonlinear optical propagation

Experimental analytical reference API: `nova::optics::nonlinear_propagation` in
`include/nova/physics/nonlinear_optical_propagation.hpp`.

## Model

For uniform continuous-wave propagation in a fixed scalar optical path,
the intensity and nonlinear phase satisfy

- dI/dz = -alpha I - beta I²
- dphi/dz = (2 pi / lambda) n2 I

Inputs in SI: vacuum wavelength (m), incident intensity (W/m²), length (m),
nonnegative linear intensity attenuation alpha (1/m), nonnegative two-photon
absorption beta (m/W), and signed Kerr coefficient n2 (m²/W).

Define L_eff = (1-exp(-alpha L))/alpha, with L_eff=L when alpha=0,
and q=beta I0 L_eff. The solution is

- I(L) = I0 exp(-alpha L)/(1+q)
- J = integral of I dz = I0 L_eff log(1+q)/q, with ratio 1 at q=0
- nonlinear phase = (2 pi / lambda) n2 J
- absorbed intensity = I0 [1-exp(-alpha L)+q]/(1+q)

The implementation uses expm1 and log1p for weak attenuation. The absorbed
quantity is an intensity difference, not power or volumetric heat generation.
No arbitrary clipping of nonfinite results is performed; invalid inputs and
detected overflow are rejected. Extreme underflow and floating-point dynamic
range remain limitations.

## Verification

Registered optics reference regression cases cover the lossless limit,
Beer-Lambert law, pure two-photon absorption, interval composition, intensity
conservation, negative Kerr response, weak absorption and invalid inputs.

During development an isolated GCC C++20 smoke build of the same numerical
implementation passed with -Wall -Wextra -Wpedantic -Wconversion -Wshadow
-Werror. An independent SciPy DOP853 integration of the intensity and
integral ODEs (seed 42, 100 logarithmically sampled parameter combinations,
rtol 2e-11, atol 1e-12) gave maximum relative discrepancy 5.368e-12.
These are limited checks, not a full repository CI run.

Run the integrated regression with:
```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release -R nova_optics_reference_tests --output-on-failure
```

## Scope

Not a nonlinear Maxwell solver, waveguide mode solver, or self-consistent
electrothermal coupling. No diffraction, dispersion, pulse dynamics,
free carriers, gain, spatial mode evolution, or temperature-dependent feedback.
Material parameters are supplied by the caller; material data are not validated
by the analytical tests.
