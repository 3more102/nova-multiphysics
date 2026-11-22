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

### Independent energy-seminorm regression

The manufactured wall test uses \(u_{exact}(x)=x(1-x)\), constant
conductivity \(k=2\), source \(f=4\), zero Dirichlet values on the
two x-plates, and insulated remaining faces. A separate test computes
\(\sum_K k\int_K |\nabla u_{exact}-\nabla u_h|^2\,dV\) **without using
the residual estimator**. Because the exact gradient is affine and the
P1 gradient is constant per tetrahedron, each integral is evaluated
exactly using the first and second tetrahedral x-coordinate moments.
It checks that one marked conforming refinement and re-solve reduces
this physical energy-seminorm error on the selected benchmark, and
that total mesh volume is preserved.

This verifies one manufactured case, not general estimator reliability
or monotonic global error reduction for arbitrary data.

The heterogeneous-material regression assigns four distinct cellwise
conductivities and verifies exact parent-to-child coefficient replication.
It also checks that every exterior face is represented exactly once by
Dirichlet classification or natural-boundary residual accounting on both
coarse and refined meshes, including implicitly insulated faces.

### Supplementary exact-norm and all-Dirichlet regressions

The additional manufactured-wall test computes continuum `L2` and `H1`
seminorm errors using the independent `compute_error_norms_3d` quadrature
on the coarse and refined meshes. It requires both to decrease for the
specified mixed-boundary quadratic wall; unlike the energy-only check,
this also evaluates the field-value norm.

A second quadratic manufactured solution `u=x^2+y^2+z^2` with `k=1`,
`f=-6` and spatially varying Dirichlet data on **every** exterior face
checks that all refined boundary-node values, including new midpoints,
are recomputed from their original spatial functions. These regression
cases are not certification or universal adaptive error guarantees.
