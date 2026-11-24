# Passive microring resonator reference
The header-only C++20 API nova/physics/ring_resonator.hpp provides the ideal single-mode add-drop ring with two reciprocal, lossless directional couplers. The all-pass topology is recovered by setting t2=1. It is not a spatial FEM solver.

Let t1,t2 be real field self-couplings, k_j=sqrt(1-t_j^2), a be the round-trip field survival, and z=exp(i phi), where phi is the round-trip phase. Both half rings have field survival sqrt(a). For unit excitation at the first input and zero excitation at the second:

- denominator D=1-a*t1*t2*z
- through amplitude H=(t1-a*t2*z)/D
- drop amplitude G=-k1*k2*sqrt(a)*exp(i phi/2)/D
- circulating power immediately after the first coupler B=k1^2/|D|^2
- independently computed absorbed power A=(1-a)*(1+a*t2^2)*B

The normalized energy identity is |H|^2+|G|^2+A=1. Critical coupling occurs at phi=0 modulo 2*pi and t1=a*t2, except for a disconnected input coupler (t1=1). The exact disconnected lossless resonance returns the no-excitation solution; arbitrary initial cavity energy is outside scope. Phase-sensitive output amplitudes depend on reference planes; the powers are 2*pi periodic.

The regression target nova_ring_resonator_tests exercises critical coupling, lossless drop, intrinsic loss, disconnected and fully cross-coupled limits, near-unity coupling, phase periodicity and a passive parameter grid. No dispersion, nonlinear feedback, backscattering or geometry-to-coupling extraction is implied.
