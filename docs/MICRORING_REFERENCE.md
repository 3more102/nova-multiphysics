# Microring resonator reference model

The C++20 header nova/physics/microring.hpp evaluates a single-input all-pass or add-drop ring with lossless reciprocal directional couplers and equal propagation halves. All amplitudes are normalized to unit incident field. t1 and t2 are real field self-transmission coefficients; a is round-trip FIELD attenuation (not power attenuation). Phase is the round-trip propagation phase in radians. Set t2=1 for an all-pass device.

Let z=exp(i phi), k_j=sqrt(1-t_j^2), and den=1-a*t1*t2*z. Then the through amplitude is (t1-a*t2*z)/den and the drop amplitude is -k1*k2*sqrt(a)*exp(i phi/2)/den. Circulating power is k1^2/|den|^2, evaluated immediately after the first coupler. Propagation power absorbed on the two segments is (1-a)*(1+a*t2^2)*circulating_power, computed independently of the output powers. The residual is 1-through_power-drop_power-absorbed_power. The amplitude reference planes retain the half-round-trip phase, so drop amplitude reverses sign after a 2*pi phase increase, while its power is periodic.

Tests cover critical coupling, lossless symmetric add-drop transfer, uncoupled input, complete round-trip extinction, 1001 phase samples compared with the analytical all-pass expression and power conservation, and invalid/singular inputs. Run ctest --test-dir build -R nova_microring_tests --output-on-failure.

Scope: a lumped single-mode steady-state optical reference, not a geometry-based FEM eigenmode or nonlinear resonator solver. No backscattering, coupler dispersion, asymmetric propagation loss, thermal detuning, or material fitting. Near-singular conditions may be rejected by the power audit. The full repository build and CI must be verified before release.
