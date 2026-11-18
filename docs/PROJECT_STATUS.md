# NOVA development checkpoint

## Active numerical checkpoint

The integration parent is [PR #17](https://github.com/3more102/nova-multiphysics/pull/17), commit `3a132bc6083adf3ac375bb48673ec60483f920d4`, verified by its [Ubuntu and Windows build](https://github.com/3more102/nova-multiphysics/actions/runs/44372586103). The feature branch `feature/transient-heat-3d` adds fixed-step P1 heat. Its Linux Release build passed 151 registered tests and all seven CTest targets. The corresponding ASan/UBSan build passed all seven with local leak detection disabled. The hosted runner is needed for the separate leak check. The next foundation milestone is stationary 2D component nullspaces for both triangle orders.

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

## Numerical and product gaps

- Stationary P1/P2 2D diffusion still requires Dirichlet data. Component-wise nullspaces and the P2 integral mean need verification before insulating-only 2D studies are accepted.
- Fixed-step transient 3D heat now passes the element-matrix, mixed-flux, temporal-order, restart and energy checks. Three cooling-cube refinements were cross-checked with independent scikit-fem assembly and SciPy sparse LU. [Verification and timings](THERMAL_3D_VERIFICATION.md)
- Scalar assembly builds triplet lists and canonicalizes them with ordered maps. Adaptive 2D thermal on this branch reassembles fixed operators during step doubling. Compare #10's operator reuse before changing that path. Serial 3D heat timings and peak resident memory are recorded through 82,944 tetrahedra. Industrial-scale qualification remains open.
- Branch-local code includes one-level 2D spatial refinement and P1/P2 transfer. P2 tetrahedra, 3D refinement, general physical groups, coupled fields and independent experimental validation remain open.
- The build links the C++ standard library and the NOVA core. No third-party numerical solver is declared by CMake. Compiler runtimes, optional product dependencies, contributor rights, redistribution terms and deliverable SBOMs still require an artifact-level audit.
- Desktop, project persistence, packaging, activation and commercial distribution remain open. See [release criteria](RELEASE_CRITERIA.md).
