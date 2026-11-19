# Conservative face-integrated heat rates (experimental)

This feature postprocesses a **steady 3D P1 tetrahedral** diffusion result, without changing the FEM solution. It addresses a distinction in the existing heat diagnostics:

- `evaluate_tagged_boundary_flux_3d` integrates the raw element gradient `-K grad(T_h)`. Across a material interface, its two one-sided normal fluxes need not coincide.
- `evaluate_diffusion_energy_balance_3d` checks global weak-form conservation (including constrained-node reactions). It does **not** assign an individual conservative rate to each tetrahedral face.
- `reconstruct_conservative_face_flux_3d` reports **one equal/opposite heat rate per interior face** and chooses heat-flow corrections to satisfy each cell's integrated source budget.

## Formulation

Let `Q_f` be signed heat flow outward from the first owner of each canonical triangular face, in W. For the second owner of an interior face, the value is exactly `-Q_f`. Let `S_c = integral_cell Q(x) dV` and let `B` be the signed cell-face incidence matrix.

1. Start with a face rate from the average of the two P1 flux vectors on an interior face. Outward Dirichlet rates start from the owner's raw P1 flux. Explicit Neumann and Robin rates are fixed to the existing degree-two face quadrature, while unassigned boundary triangles are insulated.
2. Solve a sparse graph Laplacian `B diag(A_f) B^T lambda = S - B Q_initial`, with one numerical gauge pin per connected component that has no adjustable outlet.
3. Correct each adjustable face rate by `delta Q_f = A_f B^T lambda`. This minimizes the area-weighted squared change `sum_f (delta Q_f)^2/A_f`, subject to the local balance constraints (within iterative solver tolerance).
4. Reject any closed component with incompatible total generation and prescribed outward heat flow. Do **not** silently assign heat flow to an insulated boundary.

Here `A_f` is the geometric face area. The solver uses sparse CSR conjugate gradients with Jacobi preconditioning. The result includes a per-cell signed defect and a maximum defect, which must satisfy the implementation's relative numerical acceptance threshold.

## API

```cpp
#include "nova/fem/conservative_face_flux_3d.hpp"

const auto steady = nova::fem::solve_scalar_diffusion_3d(problem);
const auto local = nova::fem::reconstruct_conservative_face_flux_3d(
    problem, steady.values);

for (const auto& face : local.faces) {
    // face.vertices, face.first_cell, face.second_cell,
    // face.area, face.outward_heat_rate
}
```

For prescribed per-cell source integrals and boundary rates from another source, use the lower-level overload accepting `TetrahedralMesh3D`, one raw `P1TetrahedralFlux3D` per tetrahedron, one integrated source rate per tetrahedron, explicitly prescribed boundary heat rates, and adjustable exterior faces.

In the high-level overload, exterior Dirichlet faces and triangles with all three explicit Dirichlet vertices are adjustable. Faces not assigned by either condition, Neumann, or Robin are insulated. Interior single-node Dirichlet constraints are not modeled as a surface outlet: a nonzero point reaction cannot be unambiguously reconstructed as a physical boundary-face flux.

## Verification and limits

Regression coverage checks affine conduction, a heated all-Dirichlet tetrahedron, a two-tetrahedron shared interface, a mixed-conductivity heated Robin/Dirichlet wall, and invalid or incompatible boundary constraints.

**Important:** This is a *conservative, face-integrated graph reconstruction*, not a conforming Raviart–Thomas/BDM H(div) vector field, a guaranteed flux error bound, or a substitute for mesh refinement. Its distribution of corrected face rates is a constrained minimum-change model, not a unique physical measurement. Exact local balance by construction does not imply a correct temperature, physical heat flux, or solver convergence. Keep global weak-residual diagnostics and analytical reference tests as separate quality gates.

Supported: stationary P1 tetrahedra with fixed meshes and existing spatial/cellwise scalar or SPD tensor conductivity, source and boundary APIs. No adaptive mesh refinement, higher-order elements, transient storage term, mechanical coupling, radiation, or certified conservative error estimator is claimed.
