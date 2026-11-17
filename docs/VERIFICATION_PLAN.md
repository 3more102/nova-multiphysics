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

- Verify CSR row offsets and sorted column indices.
- Verify duplicate triplet accumulation.
- Verify matrix-vector products against hand calculations.
- Verify dimension and index errors.
- Verify conjugate-gradient residual reduction for known SPD systems.
- Verify failure reporting for non-positive curvature.

## Gate 2: one-dimensional Poisson FEM

Benchmark:
- Domain: x in [0, 1].
- Governing equation: -u'' = 2.
- Boundary conditions: u(0) = 0 and u(1) = 0.
- Analytical solution: u(x) = x(1-x).
- Element family: linear P1 segments.
- Mesh: uniform subdivisions.

Verification:
- Compare computed nodal values with the analytical solution.
- Check the algebraic residual.
- Verify expected convergence under mesh refinement using an appropriate norm.
- Test nonzero Dirichlet values with a zero-source linear solution.
- Test invalid mesh sizes and material coefficients.

## Gate 3: two-dimensional scalar diffusion

- Unit square with a manufactured analytical solution.
- Triangular P1 and P2 elements.
- Dirichlet and Neumann boundary conditions.
- Element orientation and Jacobian checks.
- L2 and H1 error convergence studies.
- Mesh-quality sensitivity.

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
