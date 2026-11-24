# Second-harmonic generation reference

The installed header `nova/physics/shg_reference.hpp` supplies independent
lossless continuous-wave, plane-wave chi(2) benchmarks. The incident harmonic
is zero. Wavelength is the fundamental vacuum wavelength; all quantities use SI.

For angular frequency omega = 2 pi c / lambda, define
eta = 2 omega^2 d_eff^2 / (n1^2 n2 epsilon0 c^3).
Here d_eff = chi_eff^(2)/2 in the usual contracted convention, including any
effective quasi-phase-matching coefficient. Set g = L sqrt(eta I_in).

- Undepleted: I2 / I_in = g^2 sinc^2(delta_k L / 2), where sinc(x)=sin(x)/x.
  For a uniform bulk medium delta_k = 4 pi (n2-n1)/lambda.
  An effective grating wavevector may be subtracted by the caller.
- Depleted, exactly effective-phase-matched: I2 / I_in = tanh^2(g);
  I1 / I_in = sech^2(g). This is the unseeded, lossless coupled-wave solution.
  The implementation retains a stable small residual pump at high conversion.

The undepleted result deliberately permits conversion greater than one, which
indicates failure of the approximation rather than physical amplification.
The depleted reference assumes effective zero mismatch; it is not valid for
arbitrary nonzero mismatch. Both exclude walkoff, absorption, focused beams,
modal overlap, dispersion across a pulse, and full-field nonlinear Maxwell FEM.
Extreme floating-point inputs may overflow and are rejected when detected.

The CTest target nova_shg_reference_tests checks SI normalization, phase
mismatch, coherence null, intensity and length scaling, coefficient sign,
weak-conversion agreement, pump depletion, power balance and invalid inputs.
