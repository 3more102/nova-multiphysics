# Normal-incidence layer power audit

The API in `nova/physics/optical_layer_power.hpp` supplements the analytical
normal-incidence scattering solver with per-layer entering and leaving
time-averaged Poynting fluxes, normalized by incident power.

Using the exp(-i omega t) convention, E = F+B and the impedance-normalized
H = n(F-B), so normalized normal flux is Re(E conj(H))/n_incident.
The layer absorption is the difference of its two boundary fluxes.
Backward reflection is computed from the exit, and the forward amplitude
is reconstructed using only decaying propagation factors; no inverse
propagation exponential is needed for thick absorbing films.

The signed balance defect is sum(layer absorptance) minus external
absorptance (1-R-T). It is reported without clamping; consumers should
evaluate it against their numerical tolerance. The incident field is
unity. Cases with extreme field amplitudes may exceed floating-point
range and are not certified. The solver is a coherent planar,
normal-incidence, isotropic nonmagnetic reference, not a FEM volumetric
heating source or oblique-incidence model.

Reference cases cover empty, lossless, mixed absorbing/transparent,
opaque, zero-thickness and subdivided homogeneous stacks. CTest target:
`nova_optical_layer_power_tests`.
