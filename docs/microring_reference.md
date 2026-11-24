# Passive microring reference model

The C++20 header nova/physics/microring.hpp provides a scalar steady-state
two-bus add-drop ring model. For real lossless self-coupling amplitudes t1,t2,
round-trip propagation field attenuation a, and round-trip phase phi:

- through = (t1 - t2*a*exp(i*phi))/(1-t1*t2*a*exp(i*phi))
- drop = -sqrt((1-t1^2)*(1-t2^2)*a)*exp(i*phi/2)/(1-t1*t2*a*exp(i*phi))
- intrinsic absorbed power = (1-a)*(1+t2^2*a)*(1-t1^2)/|denominator|^2

All powers are normalized to unit incident power. The two half trips
have equal propagation loss and phase. There is no counter-propagating
mode or second-bus input. Setting t2=1 recovers an all-pass ring.
At resonance, t1=t2*a gives critical coupling (zero through power).
The stable denominator and direct absorption expression avoid cancellation
near high-finesse resonance. The isolated lossless ring has a nonunique
internal field but a uniquely uncoupled external response.

This is an analytical photonic circuit reference, not an electromagnetic
FEM eigenmode calculation, nonlinear solver, or dispersion extraction.
Users supply the phase and attenuation. Verification includes critical
coupling, matched lossless add-drop transfer, antiresonance, extreme
coupling limits, high finesse, periodicity and power balance.

Run with CTest: ctest --test-dir <build> -R nova_microring_reference_tests
