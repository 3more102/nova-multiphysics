# Adaptive stationary diffusion cycle on 3D P1 tetrahedra (experimental)

The public function `nova::fem::solve_adaptive_diffusion_step_3d`
connects the existing 3D stationary FEM solver, residual estimator,
Dörfler bulk marking, conforming local refinement, and P1/cell transfer.
The sequence is:

1. Solve the coarse scalar diffusion problem.
2. Evaluate element source, interior normal flux-jump, and all explicitly
   prescribed or implicitly insulated natural boundary residuals.
3. Mark the smallest prefix reaching `theta * sum(eta_K^2)`.
4. Refine all marked tetrahedra conformingly, with boundary-parent maps.
5. Interpolate coarse nodal P1 values and replicate per-cell conductivity.
6. Transfer exterior face Dirichlet/Neumann callbacks and vertex constraints.
7. Solve the refined system and recompute its residual indicators.

Use:

```cpp
#include "nova/fem/adaptive_diffusion_3d.hpp"

nova::fem::ScalarDiffusion3DProblem problem;
// Supply mesh, one positive cell_diffusivity per tetrahedron, and
// selected Dirichlet faces (with optional spatial source/Neumann flux).
const auto cycle =
    nova::fem::solve_adaptive_diffusion_step_3d(problem, 0.5);
```

Inspect `cycle.adaptation.marked_elements`,
`cycle.adaptation.refinement`, `cycle.transferred_p1_values`,
`cycle.refined_problem`, `cycle.refined_solution`, and the two
`*_indicators` objects. An identically zero estimator returns the
original problem and solution without a redundant second linear solve.

## Supported boundary and material conditions

- Straight-sided, conforming first-order tetrahedral meshes.
- Positive scalar **cellwise constant conductivity** only.
- An optional spatial source callback, or an absent (zero) source.
- Per-vertex Dirichlet constraints and facewise Dirichlet functions.
- Selected Neumann outward flux functions.
- All remaining exterior faces insulated, with their zero-flux mismatch
  included explicitly in each estimator.
- All boundary callbacks are carried to children through original-face
  ownership. On a subdivided face, the callback is evaluated at the fine
  quadrature coordinates.

Cellwise source vectors, spatial conductivity, anisotropic tensors,
Robin faces, higher-order bases, curved boundary geometry, and automatic
multicycle stopping are not implemented by this wrapper. They are rejected
rather than silently omitted. The caller may feed `refined_problem`
into a later cycle if the model remains in this supported subset.

The estimator is a **heuristic spatial error indicator**, not a certified
bound on the continuum solution. Refinement does not guarantee a lower
global discretization error or a better minimum mesh quality. Independent
benchmarks and larger-mesh performance tests remain required.

## Verification

The regressions exercise a source-loaded wall with analytic quadratic
profile, zero-error/no-refinement behavior, facewise Neumann flux
propagation, exact coarse-vertex/edge-midpoint P1 prolongation,
conservation of parent element conductivity values, and rejection of
unsupported PDE coefficients and face data. The CI matrix compiles and
runs the full suite on Ubuntu and Windows with Linux ASan/LSan/UBSan.
