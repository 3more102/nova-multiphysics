# Numerical Verification Plan

## Gate 0: build and software integrity
- Build with supported C++20 compilers.
- Run unit tests on Windows and Linux.
- Enable compiler warnings.
- Test invalid dimensions, non-finite inputs, and solver breakdown paths.
- Keep the build free of mandatory third-party numerical dependencies.

## Gate 1: sparse linear algebra
- Verify CSR duplicate triplet accumulation and matrix-vector products.
- Verify dimension and index errors.
- Verify CG solution and residual reduction for known SPD systems.
- Verify Jacobi CG rejection of non-positive or non-finite diagonal entries.
- Compare Jacobi-preconditioned and unpreconditioned CG on a scaled SPD diagonal system.
- Verify BiCGSTAB on a known nonsymmetric system with a hand-computed solution.
- Verify BiCGSTAB with Jacobi preconditioning and zero-diagonal rejection.
- Verify restarted GMRES on a known nonsymmetric system.
- Verify GMRES with Jacobi preconditioning and restart-length validation.
- Verify breakdown and iteration-limit states are not reported as convergence.
- Verify reported convergence against the recomputed true residual.

## Gate 2: one-dimensional Poisson FEM
- Domain: x in [0, 1].
- Equation: -u'' = 2.
- Boundary values: u(0) = u(1) = 0.
- Analytical solution: u(x) = x(1-x).
- Linear P1 elements on uniform meshes.
- Compare computed nodal values with the analytical solution.
- Verify a nonzero Dirichlet linear solution with zero source.

## Gate 3: two-dimensional scalar diffusion

Implemented regression cases:
- Linear harmonic patch test with nonzero Dirichlet values.
- Spatially varying positive diffusivity and source.
- Reversed triangle orientation.
- Degenerate and near-degenerate triangle rejection.
- Duplicate Dirichlet constraint rejection.
- Structured rectangular mesh topology, boundary-node and boundary-edge counts.
- Mixed Dirichlet-Neumann test for u(x,y)=x with zero volume source.
- Mixed Dirichlet-Robin test for u(x,y)=x with zero volume source.
- Rejection of Neumann/Robin conditions on interior edges and negative Robin coefficients.
- Manufactured quartic solution on 4×4, 8×8, and 16×16 meshes.
- Nodal maximum-error refinement.
- L2, H1-seminorm, and H1 error calculations.
- Exact reproduction of a linear P1 field by the error-norm calculator.
- Orientation-independent error-norm integration.

The error-norm calculator integrates the squared value and gradient errors elementwise using a Duffy-transformed 5×5 Gauss-Legendre rule.

Current limitations:
- No pure Neumann nullspace handling.
- No anisotropic diffusion tensor.
- No P2 elements.
- No automatic adaptive refinement.

Required next verification:
- Compare the error norms and observed convergence rates with an independent reference implementation.
- Add boundary flux and energy-balance checks for more general Neumann and Robin problems.
- Evaluate distorted meshes and conditioning sensitivity.

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

A benchmark is accepted only when the reference problem and assumptions are documented, automated regression tests pass, error metrics and expected convergence are stated, results are reproducible, and known limitations are recorded. Agreement with another software package is useful evidence but is not a substitute for analytical verification.
