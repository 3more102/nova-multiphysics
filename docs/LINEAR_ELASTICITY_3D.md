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

Dirichlet conditions currently prescribe **all three displacement
components at selected vertices**. At least three noncollinear fully
fixed vertices are required **per connected mesh component** to remove
all six rigid-body modes. This is a deliberate restricted API, not a
general arbitrary directional roller, symmetry, or MPC constraint model.
Duplicate constraints, invalid faces/materials/force values, and
insufficient geometric anchors are rejected.

## Public API

```cpp
#include "nova/physics/linear_elasticity_3d.hpp"

nova::physics::LinearElasticity3DProblem problem;
// Supply a validated tetrahedral mesh.
// Supply one {Young's modulus (Pa), Poisson ratio} per tetrahedron.
// Fully fix at least three noncollinear vertices on each component.
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
