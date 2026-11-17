# Numerical Verification Plan

## Purpose

NOVA must establish numerical correctness before expanding feature breadth. Each implemented feature requires automated regression tests and a documented verification case.

## Gate 0: build and software integrity

- Build with supported C++20 compilers.
- Run unit tests on Windows and Linux.
- Enable compiler warnings.
- Test invalid dimensions, non-finite inputs, and solver breakdown paths.
- Keep the build free of mandatory third-party numerical dependencies.

## Gate 1: sparse linear algebra

- Verify CSR duplicate triplet accumulation and matrix-vector products.
- Verify dimension and index errors.
- Verify conjugate-gradient solution and residual reduction for known SPD systems.
- Verify failure reporting for non-positive curvature.

## Gate 2: one-dimensional Poisson FEM

Benchmark:
- Domain: x in [0, 1].
- Governing equation: -u'' = 2.
- Boundary conditions: u(0) = 0 and u(1) = 0.
- Analytical solution: u(x) = x(1-x).
- Element family: linear P1 segments.
- Mesh: uniform subdivisions.

Automated verification currently covers:
- Nodal values against the analytical quadratic solution.
- Nonzero Dirichlet values with a zero-source linear solution.
- Invalid mesh size rejection.

Remaining:
- Explicit algebraic residual assertion.
- L2 and H1 mesh-refinement convergence rates.
- Variable-coefficient manufactured solutions.

## Gate 3: two-dimensional scalar diffusion

Implemented prototype:
- Triangular P1 mesh representation and geometric validation.
- Structured rectangular mesh generator with boundary-node and boundary-edge lists.
- Positive six-point degree-four reference-triangle quadrature.
- Jacobian-based physical basis gradients.
- Positive scalar isotropic diffusivity, including spatially varying coefficients.
- Volumetric source integration.
- Strong nodal Dirichlet elimination.
- Neumann outward-flux integration on exterior linear edges using two-point Gauss quadrature.
- Linear harmonic patch test.
- Mixed Dirichlet-Neumann test for u(x,y)=x with zero volume source.
- Manufactured solution:
  [
  u(x,y)=x(1-x)y(1-y),qquad
  f(x,y)=2[x(1-x)+y(1-y)]
  ]
  on the unit square with zero Dirichlet boundary values.
- Automated maximum-nodal-error refinement comparison on 4x4, 8x8, and 16x16 rectangular subdivisions. The test requires the error to decrease by more than a factor of 3.5 at each refinement, consistent with second-order nodal convergence for this benchmark.
- Reversed element orientation, degenerate element, invalid Neumann edge, non-positive diffusivity, and duplicate constraint tests.

Current limitations:
- No Robin boundary conditions.
- No pure Neumann nullspace handling.
- No anisotropic diffusion tensor.
- No P2 elements.
- No L2/H1 error norms or automatic adaptive refinement.

Required next verification:
- L2 and H1 error calculations and observed convergence rates.
- Boundary flux and energy-balance checks for more general Neumann problems.
- Mesh distortion and conditioning sensitivity.

## Gate 4: thermal and structural modules

Thermal:
- One-dimensional conduction.
- Layered conduction.
- Manufactured solutions.
- Energy balance.

Structural:
- Constant-strain patch test.
- Linear elasticity benchmark.
- Rigid-body mode checks.
- Reaction-force equilibrium.

## Gate 5: multiphysics coupling

- Electrothermal Joule-heating benchmark.
- Segregated versus monolithic comparison on a controlled problem.
- Coupling residual and energy-balance checks.
- Mesh and solver tolerance sensitivity.

## Acceptance policy

A benchmark is accepted only when:
1. The reference problem and assumptions are documented.
2. The implementation passes automated regression tests.
3. The error metric and expected convergence are stated.
4. Results are reproducible from committed input data.
5. Known limitations are documented.

Agreement with another software package is useful evidence but is not a substitute for analytical verification.
