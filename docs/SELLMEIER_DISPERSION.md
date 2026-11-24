# Sellmeier material dispersion reference

The header-only C++20 API `nova/physics/sellmeier.hpp` accepts a vacuum wavelength in metres and dimensionless oscillator strengths B with positive resonance wavelengths r in metres. It evaluates n² = 1 + sum B λ²/(λ²-r²) and its first two analytical wavelength derivatives.

Outputs: phase and group index (dimensionless), dn/dλ (m^-1), d²n/dλ² (m^-2), beta2 = λ³ d²n/dλ²/(2πc²) (s²/m), and D = -λ d²n/dλ²/c (s/m²). Multiply D by 10^6 to obtain ps/(nm km).

The model assumes transparent isotropic bulk material. It is not a waveguide mode solver and does not include geometric dispersion, absorption, nonlinear index, or automatic coefficient validity bands. Users must supply a fit valid for their wavelength range. Poles, nonpositive squared index and nonfinite results are rejected; near poles cancellation and conditioning may limit accuracy.

Reference: I. H. Malitson, "Interspecimen Comparison of the Refractive Index of Fused Silica," JOSA 55, 1205–1209 (1965), DOI: 10.1364/JOSA.55.001205. The included regression case uses the published fused-silica fit, not a general material database.

Run the CTest target `nova_sellmeier_tests` after building. Tests cover vacuum, fused silica, numerical derivative cross-checks, dispersion identity and invalid inputs.
