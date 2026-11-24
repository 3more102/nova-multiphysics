# Single-bus optical ring/racetrack reference

Experimental analytical model in `nova/physics/ring_resonator.hpp`.
The model uses the existing optics convention exp(-i omega t) and propagation
exp(+i phi). It is not a Maxwell FEM solver.

For real field self-coupling t and round-trip field survival a (both in [0,1]),
the lossless coupler cross amplitude is i sqrt(1-t^2). With unit input, the
through and circulating amplitudes are

    H = (t - a exp(i phi)) / (1 - t a exp(i phi))
    B = i sqrt(1-t^2) / (1 - t a exp(i phi))

The independent dissipated power fraction is (1-a^2)|B|^2. The reported
balance defect is 1 - |H|^2 - (1-a^2)|B|^2. Circulating power is an internal
normalized quantity, not an additional output port. For a waveguide with
power attenuation coefficient alpha in inverse meters, a = exp(-alpha L/2)
and phi = beta L, with round-trip length L in meters and beta in rad/m.
Critical coupling occurs at t=a<1 and phi=0 modulo 2 pi.

The model assumes one traveling mode, no reflections, a lossless coupler,
and fixed linear round-trip parameters. It does not compute waveguide
dispersion, coupling coefficients, bend loss, group index, thermal shifts,
or self-consistent Kerr effects. At t=1 the disconnected ring is assigned
zero driven circulation even at the lossless resonance. Near resonance the
denominator is evaluated without subtracting nearly equal complex numbers.

CTest target: nova_ring_reference_tests. Cases include critical coupling,
antiresonance, full coupling, lossless propagation, disconnected ring,
near-unit coupling, passivity and phase periodicity across 4961 combinations,
and invalid input rejection. No full-wave or experimental validation is claimed.
