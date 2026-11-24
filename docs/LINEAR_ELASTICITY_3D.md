# Experimental 3D P1 isotropic linear elasticity

The public API at `nova/physics/linear_elasticity_3d.hpp` introduces
the first **structural-mechanics** solve path in NOVA.

## Governing model

In the small-strain regime,

- displacement `u` is a 3-component continuous P1 tetrahedral field;
- strain `epsilon(u) = (grad(u) + grad(u)^T)/2`;
- stress `sigma = lambda tr(epsilon) I + 2 mu epsilon`;
- equilibrium `-div(sigma) = f` with applied traction `sigma*n = t`.

Each cell has a constant positive Young modulus `E` and a finite
Poisson ratio `-1 < nu < 0.5`. The solver uses the standard
isotropic Lamé parameters. Anisotropy, near-incompressibility stabilization,
plasticity, and large deformations are not claimed.

The 12×12 element matrix is assembled from the exact constant physical
shape gradients for each linear tetrahedron, with absolute Jacobian
volume. Each cell can have a constant Cartesian body-force density.
Each selected exterior triangular face can have one constant global
traction vector. Remaining exterior faces carry zero applied traction.

Dirichlet conditions can prescribe full vectors or individual global Cartesian displacement components. Each face-connected solid must suppress its six rigid-body motions. Three noncollinear fully fixed vertices satisfy this check. Compatible scalar supports can model coordinate-aligned symmetry planes. Duplicate SDK constraints, invalid material/load values, and rank-deficient supports are rejected.

## Public API

```cpp
#include "nova/physics/linear_elasticity_3d.hpp"

nova::physics::LinearElasticity3DProblem problem;
// Supply a validated tetrahedral mesh.
// Supply one {Young's modulus (Pa), Poisson ratio} per tetrahedron.
// Supply supports that remove six rigid modes on each solid.
// Optional one constant body force per element and constant face tractions.
const auto solved = nova::physics::solve_linear_elasticity_3d(problem);
```

Outputs are nodal displacement triples (metres), per-cell six-component
**tensor** strain and stress tuples ordered `xx, yy, zz, xy, yz, xz`,
and strain energy `0.5 integral sigma:epsilon dV`. The off-diagonal
strain values are **not** the engineering shear strain `gamma=2*epsilon_xy`.
The returned linear CG result includes convergence diagnostics.

## Verification cases

The added regression tests cover:

1. A nontrivial affine displacement patch on a subdivided unit cube,
   recovering all six strain/stress components and the analytical
   stored energy.
2. Rigid translation plus rotation, showing negligible strain and
   no fictitious rigid-body energy.
3. An exact uniaxial traction benchmark at `E=100 Pa, nu=0` with
   `u_x=0.1x`, `sigma_xx=10 Pa`, and `U=0.5 J` on a unit cube.
4. Reversal of tetrahedral element winding.
5. Rejection of invalid constitutive coefficients, missing materials,
   incomplete/collinear anchors, duplicate tractions, and non-finite loads.

Validation is limited to these manufactured/discrete benchmarks; no
independent general-purpose structural solver comparison, contact,
buckling, or engineering stress certification is implied.

Run the standard Windows/Linux CI matrix and Linux
ASan/LSan/UBSan jobs before marking this feature ready for merge.

## Prescribed thermal expansion (experimental thermoelastic extension)

The one-shot API accepts a **prescribed temperature change** for small-strain, cellwise isotropic thermal expansion. The prepared API also consumes nodal temperatures directly from a steady or transient heat solve.

For every tetrahedron supply both the cell_thermal_expansion_per_kelvin
array (alpha, 1/K) and the cell_temperature_change_kelvin array
(deltaT, K). Alternatively leave **both arrays empty** to retain the
original zero-thermal-strain model. Each nonempty array must have exactly
one finite value per cell. Cooling and negative thermal expansion are valid.
Mismatches, NaNs and overflowed eigenstrains/eigenstresses are rejected.

For isotropic Lamé coefficients lambda and mu:

    epsilon_th  = alpha * deltaT * I
    epsilon_el  = sym(grad(u)) - epsilon_th
    sigma       = lambda*tr(epsilon_el)*I + 2*mu*epsilon_el
    K*u         = F_body + F_traction + F_thermal
    F_thermal[a,i] = integral (3*lambda+2*mu)*alpha*deltaT*grad(N_a)[i] dV
    U           = 0.5 * integral epsilon_el : sigma dV

The cell_strain output reports **total** displacement strain; cell_stress
and strain_energy_joules use **elastic** strain, subtracting thermal
eigenstrain. Element-constant thermal loads are integrated exactly for
linear tetrahedra.

The independent analytical benchmarks use a unit cube with E=1200 Pa,
nu=0.25, alpha=0.002/K, deltaT=5 K:

1. Stress-free expansion: u=0.01(x,y,z); all stresses and elastic energy
   vanish, whereas total normal strains each equal 0.01.
2. Fully clamped heating: u=0; sigma_xx=sigma_yy=sigma_zz=-24 Pa,
   shears vanish, and elastic strain energy is 0.36 J.
3. Reject mismatched arrays, nonfinite and overflowed thermal input.
   Explicit zero-thermal-strain input must preserve the prior behavior.

The Gmsh stationary and transient heat executables can drive this mechanical solve from their computed temperatures. [Formulation, examples and measured errors](THERMOELASTIC_3D.md).
