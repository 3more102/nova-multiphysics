# Adaptive stationary diffusion cycle on 3D P1 tetrahedra (experimental)

The one-cycle function `nova::fem::solve_adaptive_diffusion_step_3d`
connects the existing 3D stationary FEM solver, residual estimator,
Dörfler bulk marking, conforming local refinement, and P1/cell transfer.
The multi-cycle driver `nova::fem::solve_adaptive_diffusion_3d` carries the
solution and indicator forward until the target or cycle limit is reached. The one-cycle sequence is:

1. Solve the coarse scalar diffusion problem.
2. Evaluate element source, interior normal flux-jump, and all explicitly
   prescribed or implicitly insulated natural boundary residuals.
3. Mark the smallest prefix reaching `theta * sum(eta_K^2)`.
4. Refine all marked tetrahedra conformingly, with boundary-parent maps.
5. Interpolate coarse nodal P1 values and transfer cellwise scalar/tensor
   conductivity and additive cell-source data.
6. Transfer Dirichlet, Neumann and Robin callbacks plus vertex constraints.
7. Solve the refined system and recompute its residual indicators.

Use:

```cpp
#include "nova/fem/adaptive_diffusion_3d.hpp"

nova::fem::ScalarDiffusion3DProblem problem;
// Supply a mesh, exactly one supported scalar/tensor material field, and
// any supported source and boundary conditions.
const auto cycle =
    nova::fem::solve_adaptive_diffusion_step_3d(problem, 0.5);
```

Inspect `cycle.adaptation.marked_elements`,
`cycle.adaptation.refinement`, `cycle.transferred_p1_values`,
`cycle.refined_problem`, `cycle.refined_solution`, and the two
`*_indicators` objects. An identically zero estimator returns the
original problem and solution without a redundant second linear solve.



## Multi-cycle adaptive driver

Use `solve_adaptive_diffusion_3d` when the mesh should be refined repeatedly
without manually wiring each refined problem into the next solve:

```cpp
nova::fem::AdaptiveDiffusionOptions3D options;
options.dorfler_theta = 0.5;
options.absolute_tolerance = 1.0e-6;
options.relative_tolerance = 1.0e-3;
options.max_cycles = 8;

const auto run = nova::fem::solve_adaptive_diffusion_3d(problem, options);
```

The stopping target is
`max(absolute_tolerance, relative_tolerance * initial_indicator)`, where
`initial_indicator` is the first coarse solve's total residual indicator.
The driver checks the initial indicator before marking. An accepted input mesh
returns its original field and one history entry with zero marked cells and
identical coarse and final sizes. Each newly refined mesh is solved and
estimated once. A failed linear solve reports the mesh stage, iteration count,
residual and CG termination message.

The result contains the final problem, solution, indicators, and a compact
per-cycle history of mesh sizes, marked elements, and indicator values.
Stop reasons distinguish a satisfied tolerance, no further refinement, and
the maximum-cycle limit. These tolerances apply to the heuristic indicator;
they do not certify the continuum solution error.

## Supported boundary and material conditions

- Straight-sided, conforming first-order tetrahedral meshes.
- Exactly one conductivity representation accepted by the base solver:
  spatial scalar, spatial symmetric positive-definite tensor, cellwise
  positive scalar, or cellwise symmetric positive-definite tensor.
- Tensor components use global Cartesian axes; all six independent symmetric
  entries, including off-diagonal terms, are retained.
- An optional spatial source callback plus an optional additive cellwise source
  vector. Cell-source values are copied from each parent tetrahedron to all
  children; the residual estimator integrates the sum of spatial and cellwise
  source densities.
- Per-vertex Dirichlet constraints and facewise Dirichlet functions.
  A boundary triangle whose three nodes are explicitly constrained (and
  which has no separately prescribed natural-boundary model) has a fully
  prescribed P1 trace. Its natural-boundary residual is excluded, and any new
  edge-midpoint DOFs on such triangles receive the linear midpoint value of
  their two prescribed parents. A face with only some constrained vertices
  remains a pointwise-constraint case and is not reclassified as fully Dirichlet.
- Selected Neumann outward-flux functions and Robin transfer/ambient functions.
- All remaining exterior faces are insulated, with their zero-flux mismatch
  included explicitly in each estimator.
- Dirichlet, Neumann and Robin callbacks are carried to children through
  original-face ownership. On a subdivided face, callbacks are evaluated at
  the fine quadrature coordinates.
- Cellwise scalar/tensor material data are copied from each parent tetrahedron
  to all of its children. Spatial material callbacks remain attached to the
  refined problem and are re-evaluated by the FEM quadrature.

Higher-order bases and curved boundary geometry remain unsupported. The
multi-cycle driver uses the same marking, boundary transfer and field estimation
without changing the estimator's heuristic status.

The estimator is a **heuristic spatial error indicator**, not a certified
bound on the continuum solution. Refinement does not guarantee a lower
global discretization error or a better minimum mesh quality. Independent
benchmarks and larger-mesh performance tests remain required.

## Verification

The regressions additionally verify nonlinear spatial face callbacks on refined nodes, exclude artificial flux residuals on fully vertex-constrained faces, and preserve the interpolated P1 trace for vertex-only boundary constraints. The regressions exercise a source-loaded wall with analytic quadratic
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


## Runnable manufactured-wall demo and VTK export

With `NOVA_BUILD_EXAMPLES=ON`, build and run
`nova_adaptive_diffusion_3d` to solve a unit cube with
`u(x)=x(1-x)`, cellwise conductivity `k=2`, source `f=4`,
Dirichlet conditions on the x-plates, and insulated other faces.
The example prints coarse/refined mesh sizes, marked-element count,
residual indicators, minimum mesh mean ratio, and an independent exact
energy-seminorm error computed using tetrahedral x-moments.

Optional ParaView export:

```sh
nova_adaptive_diffusion_3d --vtk coarse.vtk refined.vtk
```

Both files contain a first-order tetrahedral unstructured grid and the
nodal `temperature` field. A standard-library-only Python CTest verifier
runs the executable, parses both exported VTK files, checks mesh connectivity,
finite data, unit-cube volume conservation, retained Dirichlet values and
the increase in resolution. The demo fails if the independently measured
manufactured energy error does not decrease. This is a representative
verification case, not a guarantee for arbitrary adaptive refinements.


### Tensor and Robin regression coverage

The adaptive workflow now reconstructs the physical cell flux through the same
material evaluator used by the base FEM solve, then feeds those fluxes into the
residual estimator. This preserves off-diagonal anisotropic conductivity in
interior-face jumps instead of approximating a tensor by a scalar coefficient.
Regression coverage includes an affine field with a non-diagonal tensor and
its exact flux, plus a locally refined heterogeneous cell-tensor model with
both Robin and Neumann boundary faces. Parent-to-child tensor transfer and
both boundary callback transfers are checked explicitly. A further regression
combines a spatial source callback with nonuniform cellwise source densities
and verifies exact parent-to-child source transfer.

The indicator remains heuristic. For spatially varying coefficients, the
current cell residual does not include within-cell coefficient-gradient terms,
and the estimator is not a certified continuum error bound.


### Two-layer tensor wall with mixed boundaries

The unit cube has an interface at x = 0.5, a left temperature of 5.5,
a right ambient of 1 and a right transfer coefficient of 4. The left and
right conductivities have axial entries 2 and 8. Their respective first
columns are (2, 0.5, -0.25) and (8, 2, -1). The lateral faces have signed
outward fluxes from the constant vector (8, 2, -1); the volume source is zero.

For unit cross section the axial heat flux is
`(5.5 - 1) / (0.5/2 + 0.5/8 + 1/4) = 8`. The independent temperature
reference at x is `5.5 - 8*(min(x, 0.5)/2 + max(x - 0.5, 0)/8)`.
On interface-aligned grids with 2 and 4 subdivisions per axis, the exact
temperature belongs to the P1 space. At CG relative tolerances of 1e-10
and 1e-12, the regression checks every nodal temperature within 2e-8,
every cell flux component within 2e-7, and a global indicator at most 1e-7.
The benchmark uses synthetic coefficients and an exact discrete patch. The
separate quadratic manufactured wall measures discretization error under
refinement.

A source-loaded wall with relative indicator tolerance 1 checks acceptance
at equality for a nonzero initial indicator. A two-cycle regression counts
source evaluations and compares the final solution and indicator with two
explicit one-cycle calls. It verifies the work saved on their shared mesh.
The installed SDK consumer also links and executes the multicycle entry point.
