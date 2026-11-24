# Stationary 3D photothermal coupling

The one-way solver combines the existing driven complex H(curl) Maxwell solver
and continuous P1 heat conduction using the identical tetrahedral mesh.
For the exp(-i omega t) convention, passive isotropic dielectric loss produces
a time-averaged volumetric heat density of omega * eps0 * Im(eps_r) * |E|^2 / 2.
The Maxwell engine integrates that quadratic exactly on each affine cell.
Each integral Q_i (W) is divided by cell volume V_i (m^3) and added to the
thermal cell source. This preserves total absorbed power in the assembled
thermal source; it does not assert local pointwise thermal flux conservation.
Existing spatial and cell sources are retained, and input is not mutated.

Use nova::physics::Photothermal3DProblem and solve_one_way_photothermal_3d.
All mesh coordinates are in metres; temperature is in K and conductivity is
in W/(m K). Optical sources are the existing Maxwell weak PDE sources, not
incident-power-normalized ports. The result contains both solver diagnostics,
edge ordering, cell absorption, the heated thermal problem and the assembled
weak thermal energy balance. A failed solve throws rather than publishing
a misleading coupled result.

Tests target a driven constant field with independently calculated absorption,
additive thermal sources, global weak balance, zero absorption and bad cell counts.
This is a one-way stationary foundation, not a self-consistent thermo-optic,
Kerr, PML, waveguide-port, or transient optical solver. CI must pass before
treating the implementation as validated.
