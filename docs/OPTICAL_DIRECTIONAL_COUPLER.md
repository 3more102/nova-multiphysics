# Uniform optical directional coupler reference

The installed header `nova/physics/directional_coupler.hpp` provides a
closed-form reciprocal two-mode propagation benchmark. Inputs are complex
power-normalized amplitudes (sqrt(W)); length is in metres, real coupling
kappa and propagation-constant mismatch Delta beta are in rad/m, and common
power attenuation alpha is in 1/m. The average propagation phase is removed.

With d = Delta beta/2 and H = [[d,kappa],[kappa,-d]], the envelope equation is
dA/dz = i H A - alpha A/2. Since H^2 = g^2 I, g = hypot(kappa,d), the exact
propagator is exp(-alpha L/2) [cos(gL) I + i sin(gL) H/g]. The g=0 limit is
handled explicitly. A single-port input has cross-power fraction
exp(-alpha L) (kappa/g)^2 sin^2(gL). The lossless transfer ceiling is
(kappa/g)^2. For arbitrary coherent inputs, total power decays as exp(-alpha L).
Dissipated power uses expm1 for small attenuation; output power is computed
independently from the propagated amplitudes.

The model is not a geometry-to-coupling solver. It assumes uniform real
reciprocal coupling, identical attenuation in both modes, no reflection,
no radiative coupling, and no nonlinear response. The omitted common phase
must be restored externally if interfacing to an absolute-phase circuit.

The CTest target `nova_directional_coupler_tests` checks phase-matched
complete and 3 dB transfer, detuned transfer ceiling, coherent two-input
power conservation with loss, uncoupled phase, composition, reciprocity,
zero length, large attenuation and rejection of nonphysical parameters.
