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
- Verify restarted GMRES on a known nonsymmetric system with a full-size Krylov space.
- Verify GMRES convergence across multiple restart cycles using a smaller Krylov subspace.
- Verify GMRES with Jacobi preconditioning and restart-length validation.
- Verify breakdown and iteration-limit states are not reported as convergence.
- Verify reported convergence against the recomputed true residual.

## Gate 2: one-dimensional linear Poisson FEM
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
- Mixed Dirichlet-Neumann and Dirichlet-Robin tests.
- Rejection of Neumann/Robin conditions on interior edges.
- Rejection of duplicate/conflicting natural conditions on one edge.
- Rejection of negative Robin coefficients.
- Manufactured quartic solution on 4×4, 8×8, and 16×16 meshes.
- Nodal maximum-error refinement.
- L2, H1-seminorm, and H1 error calculations.
- Exact reproduction of a linear P1 field by the error-norm calculator.
- Orientation-independent error-norm integration.

Current limitations:
- No pure Neumann nullspace handling.
- No anisotropic diffusion tensor.
- No arbitrary higher-order elements beyond P2.
- No automatic adaptive refinement.

Required next verification:
- Compare error norms and observed convergence rates with an independent reference implementation.
- Add boundary flux and energy-balance checks for more general Neumann and Robin problems.
- Evaluate distorted meshes and conditioning sensitivity.

## Gate 4: nonlinear solver and nonlinear FEM
- Scalar square-root residual (F(x)=x^2-2) with analytical root (sqrt{2}).
- Verify Newton convergence from a moderate initial guess.
- Verify Armijo backtracking from a deliberately large Newton step.
- Verify non-finite trial residuals are rejected and backtracked.
- Verify both GMRES and BiCGSTAB linearized-solver paths.
- Verify singular Jacobian failure is reported without false convergence.
- Verify residual/Jacobian dimension mismatch rejection.
- Verify nonlinear P1 diffusion for (-d[(1+u^2)u']/dx=f(x)) against the manufactured solution (u=sin(pi x)).
- Verify decreasing nodal error for 8-, 16-, and 32-element meshes.
- Verify nonzero endpoint Dirichlet values.

## Gate 5: thermal and structural modules

Implemented transient thermal regression cases:
- P1 consistent mass matrix on the reference triangle, checked against exact entries.
- P2 mass symmetry, constant-field integration, and stiffness nullspace.
- Exact linear space-time temperature field with P1 and P2 Backward Euler.
- First-order temporal convergence of a quadratic-in-time manufactured P2 field.
- Second-order Crank–Nicolson temporal convergence of a cubic-in-time manufactured P2 field.
- Crank–Nicolson exact linear space-time P1/P2 solution.
- Restart at a nonzero physical start time with exact P1/P2 initial DOF values.
- Adaptive constant-field preservation and discrete thermal energy budgets.
- Adaptive Backward Euler time-dependent source and rejected-step control.
- Adaptive Crank–Nicolson tolerance tightening on a cubic-time manufactured P2 field.
- Adaptive rejection of invalid time ranges, tolerances, attempt limits and non-finite restart values.
- Global heat-content budgets under constant volumetric heating.
- Time-dependent Neumann flux integrated energy budget for Backward Euler and Crank–Nicolson.
- Time-dependent Robin convection integrated heat budget.
- Prescribed-temperature boundary reactions and discrete energy closure.
- Rejection of overlapping Dirichlet/Neumann/Robin edges, invalid Robin coefficients, and unsupported time schemes.
- Constant field with homogeneous natural heat fluxes.
- Rejection of invalid time step, negative volumetric capacity, non-finite heat source,
  interior Dirichlet edges, and conflicting corner temperature values.

Implemented nonlinear transient thermal regression cases:
- P1/P2 nonlinear temperature-dependent conductivity with an exact linear-in-space/time manufactured solution.
- Direct centered finite-difference check of the assembled P2 Jacobian applied to a test vector.
- Reversed triangle orientation with P2 nonlinear material assembly.
- Constant-temperature homogeneous-Neumann source heating.
- Linear Robin convection at the exact ambient temperature.
- Rejection of negative conductivity, non-finite constitutive temperature derivatives and incompatible initial-state DOF sizes.
- Nonconvergent Newton steps reported as failure, not a valid solution.

Required next verification:
- Independently compare nonlinear temperature-dependent conduction with a trusted reference solver.
- Add tests for large temperature gradients, distorted meshes, and strong nonlinearities.
- Verify nonlinear thermal energy and boundary flux budgets.
- Quantify discrete energy decay, adaptive error-controller sensitivity, and step efficiency for a cooling plate.
- Cross-check transient thermal benchmarks against an independent reference solver.
- Cross-check independent physical fluxes, Dirichlet boundary reactions, and energy changes.
- Add distorted-mesh verification for transient Neumann/Robin conditions.

## Legacy gate 5 cases

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

## Gate 6: multiphysics coupling
- Electrothermal Joule-heating benchmark.
- Segregated versus monolithic comparison on a controlled problem.
- Coupling residual and energy-balance checks.
- Mesh and solver tolerance sensitivity.

## Acceptance policy

A benchmark is accepted only when the reference problem and assumptions are documented, automated regression tests pass, error metrics and expected convergence are stated, results are reproducible, and known limitations are recorded. Agreement with another software package is useful evidence but is not a substitute for analytical verification.
