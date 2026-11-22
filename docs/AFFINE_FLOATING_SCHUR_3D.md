# Affine P1 floating conductors with impressed charge (experimental)

Public header: `nova/physics/affine_floating_conductor_schur_3d.hpp`.

`solve_affine_floating_conductor_schur_3d` extends the existing
source-free `solve_floating_conductor_schur_3d` without changing that API.
It supports the existing P1 linear scalar/tensor dielectric models,
optional spatial and cell free-charge densities `rho`, and prescribed
outward displacement `D.n` on exterior Neumann faces.

## Discrete equation and signs

For the FEM system `K u = F` in the interior, the terminal weak-form
reactions are **`Q_g = sum_{i in g} (K u - F)_i`**. The diffusion
assembly uses positive `+integral rho*N_i` volume loads and negative
`-integral (D.n)*N_i` exterior Neumann loads. The affine terminal law is

```text
Q = C V + Q0
C  = K_EE - K_EI K_II^{-1} K_IE
Q0 = [K u0 - F]_E
u0 = loaded P1 field with all electrode potentials grounded
C_FF V_F = Q_requested,F - Q0,F - C_FK V_K
```

`C` is extracted by the existing source-free, discrete-reaction
basis-solve API; the zero-voltage loaded solve supplies `Q0`. The
new helper integrates source and Neumann nodal loads with the same
degree-two tetrahedral and triangle quadrature rules as
`solve_scalar_diffusion_3d`. A final loaded field solve reconstructs
the potential, E and D. As before, other potential references are set
to 0 V rather than preserving their original numerical values.

The result reports the **load-corrected weak reaction** on each
conductor, the independent raw boundary-displacement estimate, and
the maximum weak-reaction mismatch on floating conductors. Do not
interpret the difference between raw boundary flux and weak
reaction as an automatic local error estimate.

## Example

```cpp
#include "nova/physics/affine_floating_conductor_schur_3d.hpp"

nova::physics::Electrostatics3DProblem problem;
// Set mesh, positive permittivity, optional free-charge sources and
// prescribed exterior displacement faces.
std::vector<nova::physics::ElectrodeBoundaryGroup3D> groups;
// Populate distinct, exterior, non-overlapping electrode face groups.
std::vector<nova::physics::FixedElectrodePotential3D> known{{10, 0.0}};
std::vector<nova::physics::FloatingElectrodeCharge3D> unknown{{20, 0.0}};
const auto solved = nova::physics::solve_affine_floating_conductor_schur_3d(
    problem, groups, known, unknown);
```

## Analytical and discrete checks

- Uniform unit cube, epsilon=2 F/m, rho=+1 C/m^3, fixed V(0)=0,
  floating Q(x=1)=0: terminal voltage 0.25 V for 1D P1 plate mesh.
  Affine zero-bias charge is -0.5 C on each plate. Final grounded
  plate reaction is -1 C, and total conductor reaction is -1 C.
- Spatial rho=1 and cell rho=1 add consistently.
- With rho=1 and prescribed D.n=+3 C/m^2 on the y=1 wall, total
  conductor reaction is +2 C on a unit cube, independently of the
  unknown floating terminal voltage.
- The source-free path agrees numerically with the previous discrete
  Schur solver, whose behavior remains unchanged.

## Scope and limitations

This is an algebraic **P1 Galerkin** terminal Schur solve, not a
monolithic saddle-point FEM method or certified equilibrated boundary
charge. Conductors must each be equipotential and cover complete,
disjoint exterior contact faces; overlapping Neumann/Dirichlet faces
are invalid. A reference is required to remove the gauge nullspace.
Nonlinear permittivity, interface/surface free-charge singularities,
Robin terms and time-dependent electrostatics are not supported.
Iterative solves and mesh refinement can affect terminal accuracy.
Run mesh-convergence studies before interpreting physical results.
