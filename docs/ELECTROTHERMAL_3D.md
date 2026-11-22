# 3D one-way electrothermal DC Joule coupling (experimental)

Public header: nova/physics/electrothermal_3d.hpp
Public API: solve_one_way_electrothermal_3d

## Equations and SI units

Electrical conduction: -div(sigma grad V) = s_I, where V is in volts,
positive-definite electrical conductivity sigma in S/m, J = -sigma grad V
is current density (A/m^2), and s_I is optional impressed current
density (A/m^3). Prescribed outward Neumann current flux is positive
out of the conductor.

Dissipation: q_J = E dot J = grad V dot sigma grad V (W/m^3),
E = -grad V. Each P1 tetrahedron has constant grad V; NOVA uses
the conductivity quadrature of its P1 solver for its averaged q_J.

Thermal conduction: -div(k grad T) = q_J + q_base, with k in W/(m K),
T in K and q_base the sum of existing thermal spatial and per-cell
heat sources. Thermal Neumann and Robin conditions are preserved.

The models require identical nodes, coordinates and element ordering.
Both scalar and SPD tensor conductivity are supported using the
underlying 3D diffusion solver. The returned thermal balance evaluates
the assembled finite-element weak-form reaction at constrained nodes,
not an independently equilibrated face flux.

## Reproduction

Build and run the nova_electrothermal_3d example. The unit-cube
reference has sigma=2 S/m, V(0)=0, V(1)=10 V, k=4 W/(m K),
T(0)=T(1)=300 K, and insulated lateral faces.

Exact continuum reference: V(x)=10x V, E_x=-10 V/m,
J_x=-20 A/m^2, q_J=200 W/m^3, total dissipated power=200 W,
and T(x)=300+25x(1-x) K. The numerical temperature is a P1
approximation, not the exact quadratic between mesh nodes.

## Limitations

Sequential one-way stationary coupling only. Conductivity is not
temperature-dependent; no nonlinear feedback, thermoelectric effects,
contact resistance, Joule heating on a second mesh, Maxwell transient
terms, thermal radiation or certified error bounds. Material and
boundary data must be supplied consistently in SI units.
