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
  A boundary triangle whose three nodes are explicitly constrained (and
  which has no separately prescribed face flux) has a fully prescribed P1
  trace. Its natural-boundary residual is excluded, and any new edge-midpoint
  DOFs on such triangles receive the linear midpoint value of their two
  prescribed parents. A face with only some constrained vertices remains a
  pointwise-constraint case and is not reclassified as fully Dirichlet.
- Selected Neumann outward flux functions.
- All remaining exterior faces insulated, with their zero-flux mismatch
  included explicitly in each estimator.
- All boundary callbacks are carried to children through original-face
  ownership. On a subdivided face, the callback is evaluated at the fine
  quadrature coordinates.

Cellwise source vectors, spatial conductivity, anisotropic tensors,
Robin faces, higher-order bases, curved boundary geometry, and certified automatic
error-tolerance stopping are not implemented by this wrapper. They are rejected
rather than silently omitted. The caller may feed `refined_problem`
into a later cycle if the model remains in this supported subset.

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

## Bounded multicycle execution

The convenience runner repeats the existing complete step on the refined
problem and stops early if Dörfler marking finds no elements:

```cpp
const auto history =
    nova::fem::run_adaptive_diffusion_cycles_3d(problem, 0.5, 3);
const auto& final_mesh = history.steps.back().refined_problem.mesh;
const auto& final_field = history.steps.back().refined_solution.values;
```

`max_cycles` must be at least one. A zero-error indicator stops after a
single non-refining step; otherwise the runner performs at most the specified
number of steps. `stopped_without_marking` reports that early-stop reason.
Every step retains its coarse/refined solutions and indicators so numerical
changes can be inspected. Each coarse system is currently solved again, even
when its PDE is identical to the preceding refined system; solver-state reuse
is future work. The runner has no automatic error-tolerance, cell-budget,
or convergence guarantee: choose a small cycle cap for large 3D meshes.
