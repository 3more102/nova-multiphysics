# Coherent planar optical two-port reference

`nova/physics/coherent_two_port.hpp` extends the isotropic oblique stack reference
to simultaneous coherent illumination from two propagating, lossless exterior
ports. Both excitations share the same transverse wavevector and TE or TM
polarization. The right angle follows Snell's law; the reversed layer sequence
provides right-incidence amplitudes. Port planes lie at the outer stack faces.

Tangential electric amplitudes are normalized by the square root of real exterior
admittance (TE: n cos(theta), TM: n/cos(theta)). The outgoing power waves obey
b = S a. Reciprocity requires S12 = S21. Passivity requires I - S* S to be
positive semidefinite, where S* is the conjugate transpose. Its two eigenvalues
bound the absorbed fraction for any unit-power coherent excitation. The
implementation audits both conditions and rejects unresolved numerical states.
Only roundoff-scale eigenvalue excursions are clipped.

`coherent_power` returns incident, outgoing and signed residual absorbed power.
The residual is not an independent volume integral. A maximum near one suggests
a coherent perfect absorption condition, not proof of a root or an optimized
device. No evanescent exterior ports, anisotropy, gain, finite beams or
incoherent phase averaging are supported. Near-grazing/critical or ill-conditioned
states may be rejected.

CTest target: `nova_coherent_two_port_tests`. Tests include lossless unequal
exterior media, TE/TM reciprocity, asymmetric passive stacks, one-sided parity,
phase sweeps, independent even/odd Fabry-Perot eigenchannels and invalid ports.
