# Transient 3D heat verification

This dossier covers scalar and tensor P1 heat diffusion on affine tetrahedra with spatially variable, time-independent volume materials. The cooling-cube comparison uses fixed steps. The tensor study also checks adaptive trials and an independently integrated terminal state. Material calibration and experimental validation remain open.

## Element and boundary evidence

For a tetrahedron with coordinate intercepts 2, 3 and 4, volume is 4 and the mean of $k=1+x+2y+3z$ is 6. The test checks all 16 entries of both the exact consistent mass and gradient stiffness matrices, then reverses cell orientation.

On a distorted unit cube the manufactured field is $T=2+x+2y-z+t$, with $c=1+x$, $k=2+x$ and $Q=x$. Dirichlet data apply at $x=0$. At $x=1$ one half of the face uses Neumann data and the other uses Robin exchange with $h=2+t$. The other four faces use the analytical outward flux. Heat content is $55/12+3t/2$. The expected input per unit time is $1/2$ from the volume source, $3/2$ from Neumann faces, $3/2$ from Robin faces and $-2$ from the Dirichlet boundary. Both schemes recover the analytical temperatures and all four exchanges to the regression tolerances.

For a reference tetrahedron with its $x=0$ face held at zero and the remaining faces insulated, the free P1 degree of freedom satisfies $u'=-10u$. For unit initial amplitude, each scheme is checked against its discrete amplification factor and against $e^{-2}$ at time 0.2. At 10, 20 and 40 steps the successive error orders are 0.976 and 0.988 for Backward Euler, and 2.003 and 2.001 for Crank–Nicolson.

The GCC 14.2 Linux Release build passed 151 registered tests and seven CTest targets. ASan and UBSan passed all seven targets locally with leak detection disabled. With local leak detection enabled the sanitizer runtime terminated at exit because its thread-inspection facility was unavailable. Clang-Tidy 23.1.4 completed six changed numerical translation units with its core/C++ analyzer checks, use-after-move, dangling-handle and range-copy checks enabled and no diagnostics. The hosted Linux runner subsequently passed with leak detection enabled. The same run also passed the independent solver, static checks, Ubuntu Release and Windows Release. [Exact thermal CI run](https://github.com/3more102/nova-multiphysics/actions/runs/44375598672)

## Independent spatial comparison

The unit cooling cube has $c=k=1$, zero Dirichlet temperature on all faces, zero source and initial temperature $\sin(\pi x)\sin(\pi y)\sin(\pi z)$. The analytical amplitude at time $t$ is $e^{-3\pi^2 t}$. Both solvers use 20 Crank–Nicolson steps to $t=0.02$. The reference script constructs Kuhn tetrahedra independently, assembles P1 operators with scikit-fem 12.0.2 and factors the interior system with SciPy 1.17.0 sparse LU. It compares every final vertex and evaluates the continuum $L^2$ error by degree-five tetrahedral quadrature. NumPy version was 2.3.5. [Reference backend](https://scikit-fem.readthedocs.io/en/stable/)

| Subdivisions per axis | Tetrahedra | DOFs | Maximum absolute difference from reference FEM | NOVA continuum L2 error | Refinement order |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 162 | 64 | 6.11e-16 | 8.807897e-2 | — |
| 6 | 1,296 | 343 | 3.34e-12 | 2.749026e-2 | 1.6799 |
| 12 | 10,368 | 2,197 | 1.40e-12 | 7.299464e-3 | 1.9131 |

The regression limit for the backend discrepancy is 1e-10. Each refinement order must exceed 1.6 on this benchmark. The backend comparison covers this fixed-coefficient, prescribed-temperature case. Mixed-face verification uses the separate analytical problem above.

## Serial solve measurements

Linux x86-64, Debian 13, AMD EPYC 9V74, GCC 14.2, CMake Release. One process pinned to a logical CPU, one warm-up followed by five runs per size. The table gives medians. Each run uses ten Crank–Nicolson steps, Jacobi CG (relative 1e-12, absolute 1e-15) and retains all step solutions. The timed interval includes validation, assembly and stepping. It excludes structured mesh generation and report printing. Peak resident memory covers the entire example process.

| Subdivisions per axis | Tetrahedra | DOFs | Total CG iterations | Solve time (s) | Peak RSS (MiB) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 384 | 125 | 60 | 0.00109 | 2.0 |
| 8 | 3,072 | 729 | 120 | 0.01035 | 5.0 |
| 16 | 24,576 | 4,913 | 180 | 0.09069 | 27.9 |
| 24 | 82,944 | 15,625 | 260 | 0.35572 | 88.5 |

These timings describe a single serial configuration. Hardware comparisons and industrial-scale qualification require further measurements.

## Reproduction

From the repository root on Linux (use the Release executable path on Windows):

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/nova_transient_heat_3d 24 10
python -m pip install "numpy==2.3.5" "scipy==1.17.0" "scikit-fem==12.0.2"
python verification/compare_transient_heat_3d.py ./build/nova_transient_heat_3d
~~~

The last two commands are for an optional verification environment. The C++ library and examples build with the standard library alone. For Linux memory reproduction, run the example under GNU time with format `%M` for peak KiB. The example's `solve_seconds` reports the solve-only interval.

## Adaptive tensor heat

The unit cube has unit capacity, zero prescribed temperature on its exterior, source $(1+t)(1+z)$ and initial field $64x(1-x)y(1-y)z(1-z)(1+x/2)$. Conductivity is $(1+x/2)$ times a Cartesian tensor with diagonal $(2,1.5,1)$ and off-diagonal $(K_{xy},K_{xz},K_{yz})=(0.4,-0.25,0.3)$. NOVA integrates to $t=0.02$ with Jacobi CG. The adaptive relative tolerance is zero.

The reference script builds its own mesh, stiffness, mass and source using scikit-fem. SciPy sparse direct solves repeat every full and half trial from the last independently accepted state. The script checks each normalized error, acceptance decision, physical time and final nodal value. A second reference integrates the same spatial system with SciPy Radau on its own time grid at relative tolerance $10^{-11}$ and absolute tolerance $10^{-13}$. This comparison measures time-integration error at fixed spatial resolution. [SciPy Radau](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.Radau.html)

| Cells per axis | Scheme | Absolute step tolerance | Accepted / rejected | Maximum nodal difference, same steps | Maximum nodal difference, Radau |
| ---: | --- | ---: | ---: | ---: | ---: |
| 3 | Backward Euler | 1e-4 | 36 / 3 | 6.38e-16 | 3.63e-3 |
| 3 | Backward Euler | 1e-6 | 358 / 4 | 1.61e-15 | 3.69e-4 |
| 3 | Crank–Nicolson | 1e-4 | 8 / 3 | 3.89e-16 | 3.93e-4 |
| 3 | Crank–Nicolson | 1e-6 | 32 / 3 | 3.05e-16 | 2.08e-5 |
| 4 | Backward Euler | 1e-4 | 36 / 3 | 1.28e-13 | 4.64e-3 |
| 4 | Backward Euler | 1e-6 | 355 / 4 | 6.14e-12 | 4.64e-4 |
| 4 | Crank–Nicolson | 1e-4 | 8 / 3 | 1.78e-14 | 3.70e-4 |
| 4 | Crank–Nicolson | 1e-6 | 32 / 3 | 6.76e-14 | 2.16e-5 |

The repeated-step comparison accepts at most 2e-11 absolute nodal error. Tightening the tolerance by 100 must reduce the Radau discrepancy by at least four on both meshes. The tight-run terminal-error limits are 1e-3 for Backward Euler and 5e-5 for Crank–Nicolson. The mixed-face analytical case and a stationary two-layer flux-continuity case in `tests/test_anisotropic_3d.cpp` exercise boundary and interface terms separately.

### Adaptive serial measurements

Linux x86-64, AMD EPYC 9V74, GCC 14.2 Release. Each Crank–Nicolson run used one pinned logical CPU and absolute tolerance 1e-4. The medians cover three new processes after one warm-up at each size. Wall time includes launch, mesh generation, assembly, all accepted and rejected trials, and exit. Peak RSS is per process.

| Tetrahedra | DOFs | Accepted / rejected | Wall time (s) | Peak RSS (MiB) | Maximum step energy imbalance |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 384 | 125 | 8 / 3 | 0.0051 | 2.1 | 4.27e-16 |
| 3,072 | 729 | 8 / 4 | 0.0272 | 5.1 | 3.04e-15 |
| 24,576 | 4,913 | 9 / 5 | 0.2406 | 28.8 | 2.64e-15 |
| 82,944 | 15,625 | 9 / 5 | 1.1196 | 89.6 | 1.02e-14 |

~~~sh
python verification/compare_adaptive_heat_3d.py ./build/nova_anisotropic_heat_3d
./build/nova_anisotropic_heat_3d cn 24 1e-4
~~~
