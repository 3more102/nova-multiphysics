# Restricted Rayleigh sphere reference

The header-only C++20 API `nova/physics/rayleigh_sphere.hpp` is an electric-dipole reference for passive, homogeneous, nonmagnetic spheres in a lossless nonmagnetic host. Inputs are vacuum wavelength, radius (SI metres), real positive host index, and complex permittivity contrast e = epsilon_particle / epsilon_host. The time convention is exp(-i omega t).

## Model

Let x = 2*pi*n_host*a/lambda0 and g = 2*x^3/3. Define D = e+2-i*g*(e-1) and F = (e-1)/D. The radiation-reaction correction ensures the dipole optical theorem, including for lossless particles.

- Q_sca = (8/3)*x^4*|F|^2.
- Q_abs = 12*x*Im(e)/|D|^2.
- Q_ext = Q_sca + Q_abs = 4*x*Im(F).
- Cross sections in square metres equal efficiencies multiplied by pi*a^2.

The material absorption expression avoids subtracting two nearly equal terms. For a lossless sphere Q_abs is exactly zero. For small nonresonant particles scattering cross section approaches radius^6 scaling.

## Scope

The routine rejects x > 0.1, |sqrt(e)|*x > 0.1, gain, invalid inputs, and |e+2| <= 10*x^2*max(1,|e|). The last condition is a conservative exclusion of the electrostatic pole; it is not a rigorous error estimate. The model omits dynamic depolarization, higher multipoles, nonspherical geometry, substrate coupling, near fields and full Mie scattering. Do not treat it as a predictive model at a plasmonic resonance. No material database is implied.

## Verification target

`nova_rayleigh_sphere_reference_tests` exercises index matching, host wavenumber, the closed-form corrected dielectric dipole, lossless absorption, radius scaling, passive metal absorption, independent optical-theorem identity, SI units and invalid/out-of-scope inputs.

After configuring and building with NOVA's default preset, run:

```sh
ctest --preset default -R nova_rayleigh_sphere_reference_tests --output-on-failure
```

This is an analytical benchmark for future 3D Maxwell scattering work, not an implementation or validation of that solver. A full Mie series and independently converged cross-solver comparisons remain necessary.
