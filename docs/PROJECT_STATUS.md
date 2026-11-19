# NOVA development checkpoint

## 19 November 2026 · linear-time exact CSR factor-cache matching (pending CI)

The IC(0) cache now compares canonical CSR row offsets, column indices and
coefficient arrays directly, replacing repeated binary searches for individual
matrix elements. Equality remains exact; coefficient changes, dimension changes
and sparsity-pattern changes still force refactorization. The comparison costs
O(rows + nnz), with no fingerprint collisions or tolerance-based stale reuse.
A regression compares equivalent reordered triplets against changed numeric,
pattern and dimension cases. Verification on the exact branch commit is required
before this optimization is marked validated.

## 19 November 2026 · cross-trial adaptive 3D IC(0) reuse (pending CI)

This increment follows draft [PR #23](https://github.com/3more102/nova-multiphysics/pull/23).
The adaptive P1 tetrahedral step-doubling controller now owns separate IC(0)
factor caches for its full-step and half-step trial systems. The factor is
reused across distinct adaptive trials only when the *complete constrained CSR
matrix* is identical; changing dt or Robin transfer coefficients forces a new
factor. Each cache is confined to one simulation and uses the existing
true-residual convergence check. The result exposes separate factor-build
counts for both slots; Jacobi and unpreconditioned solves report zero.

Regression additions check steady full/half-step reuse for Backward Euler and
Crank–Nicolson, invalidation when the last macro-step has a different dt, and
invalidation with time-varying Robin exchange. The repository branch must pass
the complete Windows/Linux, sanitizer and independent-FEM CI jobs before this
checkpoint can be labeled verified. This does not merge or supersede PR #23
or the upstream stack.

## 19 November 2026 · cached IC(0) for repeated thermal matrices

The CG API now includes `CachedConjugateGradient`, an instance-owned IC(0)
factor cache. The cache compares every stored nonzero against the supplied
canonical CSR matrix using exact coefficients, dimensions and sparsity count.
It rebuilds the incomplete factor if any entry or pattern changes, so updates
to timestep, conductivity/capacity operators, constrained rows or Robin
transfer coefficients cannot silently reuse a stale factor. The same factor
can solve successive right-hand sides with different prescribed temperatures
and heat sources. Non-IC(0) solver choices retain the uncached code path.

The shared fixed-step thermal integrator owns one cache for its integration
lifetime. The no-Robin constrained matrix already had assembly reuse; now
its IC(0) factor is also reusable. When Robin coefficients vary, the full
matrix is compared and factorization is refreshed as necessary. Adaptive
3D step-doubling receives reuse within each two-half-step integration, but
the cache is not yet shared across separate adaptive trials. Solver
termination still verifies the independently recomputed true residual.

New tests cover repeated RHS, changed numerical coefficients, changed sparsity,
switching between Jacobi and IC(0), and agreement on a 3D heat problem with
time-dependent Robin coefficients. Verification on this branch must pass CI
before the increment is treated as validated.


## 19 November 2026 · zero-fill incomplete Cholesky

The solver increment builds on the verified tensor branch at `92960e071255`. Select `PreconditionerType::IncompleteCholesky0` in the existing CG options for stationary or transient problems. The public factor also supports repeated triangular solves. Seven additional regressions cover a four-cycle with known fill, exact one-iteration CG and its true residual, scaling, nonsymmetric input, breakdown on a verified SPD matrix, an anisotropic stationary patch and adaptive heat history. GCC 14.2 Release passed all 177 unit tests and all nine CTest targets. [Parent Windows/Linux, sanitizer and reference CI](https://github.com/3more102/nova-multiphysics/actions/runs/44478548702)

The stationary cubic solution uses unit scalar diffusivity, prescribed zero temperatures on both x-faces, insulated lateral faces and CG relative tolerance 1e-11. The example reports the Euclidean true residual and maximum nodal error. Wall time covers the complete process, including mesh generation, assembly and factor setup. Medians use five alternating runs per method after one warm-up per method on one pinned AMD EPYC 9V74 logical CPU.

| Tetrahedra | DOFs | Iterations, Jacobi / IC(0) | True residual, Jacobi / IC(0) | Max nodal error, both | Wall time (s), Jacobi / IC(0) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3,072 | 729 | 50 / 21 | 7.74e-13 / 3.28e-13 | 2.55202e-3 | 0.0071 / 0.0071 |
| 24,576 | 4,913 | 98 / 39 | 3.80e-13 / 2.76e-13 | 6.53796e-4 | 0.0494 / 0.0503 |
| 82,944 | 15,625 | 142 / 56 | 1.96e-13 / 1.73e-13 | 2.91934e-4 | 0.2143 / 0.2211 |

Reproduce the last row with `nova_diffusion_3d jacobi 24` and `nova_diffusion_3d ic0 24`. IC(0) reduced iterations on every mesh. Its total time on the largest case was about 3% higher. The next solver milestone is reuse of a prepared factor across identical heat-step matrices, with explicit invalidation on step size, material, constraint or Robin-matrix changes.

## 19 November 2026 · tensor diffusion and adaptive thermal verification

The pre-change audit found 19 open PRs. All current base commits were ancestors of their PR heads. The default head was `0d56a5efe268`, and the full tetrahedral stack ended at #20 (`b68cbd326a4b`). Its clean Linux Release build passed 163 unit tests and eight CTest targets. The stationary 3D example used 3,072 tetrahedra, 50 CG iterations and had a true residual of 7.74e-13 with maximum nodal error 0.00255202. The dependency table below identifies the six open sibling PRs outside the 3D path.

[PR #21](https://github.com/3more102/nova-multiphysics/pull/21) targets #20. The shared P1 kernel now accepts Cartesian SPD conductivity for steady, fixed-step and adaptive solves. Two scalar range regressions reproduced on #20. A conductivity of 1e-90 on a right tetrahedron of side 1e-90 gave a zero thermal stiffness entry. A stationary patch with conductivity 1e200 on the same mesh raised overflow. The corrected kernel returns the analytical 1.66667e-181 and the patch value 1.0. The implementation passed seven new tests and all four hosted jobs. [Implementation CI](https://github.com/3more102/nova-multiphysics/actions/runs/44477873595)

The adaptive comparison covers eight combinations of mesh, scheme and tolerance. The Linux suite including the new example passed 170 unit tests and nine CTest targets. The largest same-step difference from independent sparse direct solves was 6.15e-12 or less. The separate time-adaptive Radau comparison and serial timings through 82,944 tetrahedra are recorded in the [verification dossier](THERMAL_3D_VERIFICATION.md#adaptive-tensor-heat). Remaining work includes stronger sparse preconditioning, higher-order tetrahedra, material calibration and independent physical validation.

## Numerical checkpoint through 18 November 2026

| Checkpoint | Tested scope | Evidence |
| --- | --- | --- |
| [#17](https://github.com/3more102/nova-multiphysics/pull/17), `3a132bc6083a` | 3D diffusion integration, mesh quality and diagnostics | 141 registered tests and six CTest targets. [Ubuntu and Windows](https://github.com/3more102/nova-multiphysics/actions/runs/44372586103) |
| [#18](https://github.com/3more102/nova-multiphysics/pull/18), `7a055a401219` | Transient P1 tetrahedral heat and continuous numerical verification | 151 registered tests and seven CTest targets. [Ubuntu, Windows, ASan with leak detection, UBSan, independent FEM and static analysis](https://github.com/3more102/nova-multiphysics/actions/runs/44375598672) |
| `feature/component-neumann-2d` on #18 | Stationary component-wise 2D P1/P2 nullspaces and triangle topology | 158 registered tests and seven CTest targets passed in the local GCC 14.2 Release build. The independent 3D reference comparison still passes. |
| [#20](https://github.com/3more102/nova-multiphysics/pull/20), verified source `69be028397a6` | Adaptive P1 tetrahedral Backward Euler/Crank–Nicolson step doubling, static operator reuse, exact restarts and energy diagnostics | Five new regression cases; all eight CTest targets pass on Ubuntu, Windows and ASan/UBSan. Static analysis and the independent fixed-step FEM comparison pass. [Exact-head CI run 44379306942](https://github.com/3more102/nova-multiphysics/actions/runs/44379306942). |

The adaptive 3D milestone is implemented and CI-verified on `feature/adaptive-heat-3d`: full-step/two-half-step error control, accepted/rejected trials, component heat preservation, exact restart vectors and one-time static operator assembly. PR #20 remains a draft stacked on #19; no upstream changes were merged. Independent adaptive cross-solver validation and experimental validation remain open. Next thermal priorities are nonlinear material properties and higher-order tetrahedral spaces.

## Baseline audit · 18 November 2026

CMake declares version 0.8.0. The default branch at audit was `0d56a5efe268`. All 15 open PRs had their current base commit as an ancestor of their head. The 3D diffusion head `1a691dbeb091` was 85 commits ahead of main with zero commits missing from main. On an independent Linux Release build with GCC 14.2, this head passed 138 registered tests and all six CTest targets. The diffusion example used 3,072 tetrahedra, 729 DOFs and 50 CG iterations. It reported a true residual of 7.73738e-13 and a maximum nodal error of 0.00255202.

The table records the upstream PRs **before** the integration checkpoint. Each linked run is a successful push run on the specified branch and exact commit. Stacked PRs do not trigger the workflow's pull-request event unless they target main.

| PR | Base | Head at audit | Scope | Exact-head CI |
| --- | --- | --- | --- | --- |
| [#1](https://github.com/3more102/nova-multiphysics/pull/1) | `main` | `f684147a` | Triangle H1 P1/P2 | [run 44255250308](https://github.com/3more102/nova-multiphysics/actions/runs/44255250308) |
| [#2](https://github.com/3more102/nova-multiphysics/pull/2) | #1 | `fdc3c17f` | 2D consistent mass and Backward Euler | [run 44256915422](https://github.com/3more102/nova-multiphysics/actions/runs/44256915422) |
| [#3](https://github.com/3more102/nova-multiphysics/pull/3) | #2 | `9a9b391b` | Crank–Nicolson and thermal boundaries | [run 44258282288](https://github.com/3more102/nova-multiphysics/actions/runs/44258282288) |
| [#4](https://github.com/3more102/nova-multiphysics/pull/4) | #3 | `d81b7f96` | Adaptive time stepping | [run 44259273586](https://github.com/3more102/nova-multiphysics/actions/runs/44259273586) |
| [#5](https://github.com/3more102/nova-multiphysics/pull/5) | #4 | `b1353dc3` | VTK field export | [run 44260278170](https://github.com/3more102/nova-multiphysics/actions/runs/44260278170) |
| [#6](https://github.com/3more102/nova-multiphysics/pull/6) | #4 | `177906d6` | Nonlinear Picard heat | [run 44260536808](https://github.com/3more102/nova-multiphysics/actions/runs/44260536808) |
| [#7](https://github.com/3more102/nova-multiphysics/pull/7) | #4 | `4685edad` | Nonlinear Newton heat | [run 44264595218](https://github.com/3more102/nova-multiphysics/actions/runs/44264595218) |
| [#8](https://github.com/3more102/nova-multiphysics/pull/8) | #4 | `5fbe53c8` | Triangle quality and uniform refinement | [run 44261254989](https://github.com/3more102/nova-multiphysics/actions/runs/44261254989) |
| [#9](https://github.com/3more102/nova-multiphysics/pull/9) | #8 | `aa72a80b` | P1/P2 transfer and remeshing | [run 44263979013](https://github.com/3more102/nova-multiphysics/actions/runs/44263979013) |
| [#10](https://github.com/3more102/nova-multiphysics/pull/10) | #4 | `4e8d21ec` | Adaptive operator cache | [run 44263089382](https://github.com/3more102/nova-multiphysics/actions/runs/44263089382) |
| [#12](https://github.com/3more102/nova-multiphysics/pull/12) | #9 | `0dbc94e8` | Red-green adaptivity | [run 44264847809](https://github.com/3more102/nova-multiphysics/actions/runs/44264847809) |
| [#13](https://github.com/3more102/nova-multiphysics/pull/13) | `feature/local-adaptive-mesh-refinement` | `26185630` | P1 residual estimator on alternate local refinement | [run 44264603272](https://github.com/3more102/nova-multiphysics/actions/runs/44264603272) |
| [#14](https://github.com/3more102/nova-multiphysics/pull/14) | #12 | `99593222` | P1/P2 residual estimator on red-green refinement | [run 44265026039](https://github.com/3more102/nova-multiphysics/actions/runs/44265026039) |
| [#15](https://github.com/3more102/nova-multiphysics/pull/15) | #12 | `3f2f354b` | Tetrahedral geometry and H1 | [run 44265744204](https://github.com/3more102/nova-multiphysics/actions/runs/44265744204) |
| [#16](https://github.com/3more102/nova-multiphysics/pull/16) | #15 | `1a691dbe` | P1 stationary 3D diffusion | [run 44368643231](https://github.com/3more102/nova-multiphysics/actions/runs/44368643231) |

## Integration and review

The 3D path is main → #1 → #2 → #3 → #4 → #8 → #9 → #12 → #15 → #16. A change can be tested on the complete ancestry by branching from #16. Maintainers must recheck the actual bases and CI before integrating into main. In particular, squashing an upstream PR changes the ancestry expected by its children.

PRs #5, #6, #7, #10, #13 and #14 are outside the 3D path. #13 and #14 overlap in residual estimation and use different refinement implementations. #13 targets a branch without an open parent PR. #11 was closed as superseded. The unproposed `feature/p1-p2-residual-estimator` and the closed #11 branch are outside the integration baseline.

| Review finding | Disposition at baseline |
| --- | --- |
| #15 repeated whole-mesh validation per element gradient | Already corrected in `3f2f354b`. The local gradient evaluates the element's four vertices. The review thread remained open at audit. |
| #15 scale dependence of dihedral angles | Reproduced on the right tetrahedron at scale 1e-90. The integration branch uses dimensionless edges and unit normals. |
| #3 Backward Euler evaluates initial loads with zero weight | Reproduced by inspection on #16 and corrected in the shared fixed-step integrator on the thermal branch. Both dimensional APIs have regression coverage. |
| #6 first Picard material evaluation before prescribed initial temperatures | Open review on a sibling branch. Verify start-time projection there before integration. |

Three new regressions failed on the original diffusion head and passed after the corrections. They cover scale/orientation invariance, folded faces and unused vertices, and correspondence between the reported CG residual and the actual gauged iterate. Linux Release then passed 141 registered tests and six CTest targets. Mesh validation checks adjacent cells across shared faces. Global intersection detection for geometrically overlapping, topologically unrelated tetrahedra remains open.


## Current open-PR dependency audit · 18 November 2026

The repository has 19 open PRs. At this audit, all 19 current base commits are ancestors of their corresponding heads (behind_by = 0), and GitHub reports each PR as mergeable. This is an ancestry/conflict check, not merge approval. The default branch remains main at 0d56a5efe268.

| Dependency | Open PRs |
| --- | --- |
| main → #1 → #2 → #3 → #4 | Core H1, transient heat, Crank–Nicolson and 2D adaptive stepping |
| #4 → #5, #6, #7, #8, #10 | VTK export, Picard thermal, Newton thermal, triangle refinement, and the separate 2D operator cache |
| #8 → #9 → #12 | P1/P2 transfer, transient remeshing and red-green adaptive refinement |
| #12 → #14 | P1/P2 residual-based diffusion estimator |
| #12 → #15 → #16 → #17 → #18 → #19 → #20 | 3D mesh/H1, stationary diffusion, integration corrections, fixed-step 3D heat, 2D component nullspaces, and adaptive 3D heat |
| feature/local-adaptive-mesh-refinement → #13 | Separate P1 residual estimator; no open parent PR |

PRs #13 and #14 cover overlapping residual-estimator functionality on different refinement branches. Review their assumptions and implementation differences before integrating both. PR #10 is the 2D adaptive operator-cache branch; PR #20 independently reuses static 3D thermal operators and does not depend on #10.

The verified 3D integration chain is #15 → #16 → #17 → #18 → #19 → #20, with #15 descending from #12. Merge or retarget in dependency order. Squashing an upstream PR changes the ancestry expected by its children. No PR in this stack was merged during this audit.

## Numerical and product gaps

- Stationary P1/P2 triangles and P1 tetrahedra now share per-component compatibility and normalization. Tests cover a floating component next to Robin- and Dirichlet-anchored components. Pinning is a simple gauge choice. Conditioning and preconditioning on large floating systems still need study.
- Fixed-step transient 3D heat now passes the element-matrix, mixed-flux, temporal-order, restart and energy checks. Three cooling-cube refinements were cross-checked with independent scikit-fem assembly and SciPy sparse LU. [Verification and timings](THERMAL_3D_VERIFICATION.md)
- Scalar assembly builds triplet lists and canonicalizes them with ordered maps. Adaptive 2D thermal on this branch reassembles fixed operators during step doubling. Compare #10's operator reuse before changing that path. Serial 3D heat timings and peak resident memory are recorded through 82,944 tetrahedra. Industrial-scale qualification remains open.
- Branch-local code includes one-level 2D spatial refinement and P1/P2 transfer. P2 tetrahedra, 3D refinement, general physical groups, coupled fields and independent experimental validation remain open.
- The build links the C++ standard library and the NOVA core. No third-party numerical solver is declared by CMake. Compiler runtimes, optional product dependencies, contributor rights, redistribution terms and deliverable SBOMs still require an artifact-level audit.
- Desktop, project persistence, packaging, activation and commercial distribution remain open. See [release criteria](RELEASE_CRITERIA.md).
