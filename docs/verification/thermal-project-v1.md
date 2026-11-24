# Stationary thermal box verification, schema 1

| Run identity | Value |
| --- | --- |
| Solver implementation | `0050428fa149848a552b85c66f7c74985ef2d4b6` |
| Starting main | `0a63e2460bef2506bc7ea4137d61abbbd464ce40` |
| Merged main | `6a86402ede6238e810b50c7e96fbba25da844bcb` (same Git tree) |
| Local build | Debian 13 x86-64, GCC 14.2, CMake 3.31.6, Release; examples, tests and Python verification enabled |
| CTest | 41 passed, 0 failed, 0 skipped; includes 396 cases inside the main C++ test executable |
| Installed-package consumer | 4 passed, 0 failed, 0 skipped |
| Hosted checks | Ubuntu Release, Windows Release, Linux address/leak/undefined-behavior checks passed. [PR workflow](https://github.com/3more102/nova-multiphysics/actions/runs/45443673158) |

## Exact wall with convection

| Quantity | Analytical reference | Absolute acceptance | Measured |
| --- | ---: | ---: | ---: |
| Cold face | 325 K | 2e-8 K | 324.99999999995015 K |
| Convective outward rate | 75 W | 2e-8 W | 75.00000000000014 W |
| Fixed-temperature outward rate | −75 W | 2e-8 W | −74.99999999999885 W |
| Maximum nodal error | 0 K | 2e-8 K | 4.98517e-11 K |
| Maximum cell vector-flux error | 0 W/m² | 1e-6 W/m² | 9.19457e-10 W/m² |
| Signed heat balance defect | 0 W | 2e-8 W | −1.29873e-12 W |

## Independently integrated heated-wall solution

A uniform 240 W/m³ source generates 30 W. The specified outward flux removes 10 W at the right boundary. The fixed-temperature reaction removes 20 W. The checker evaluates the independently derived quadratic temperature and linear Fourier heat flux at cell centroids. Each refinement must exceed orders 1.6 (temperature) and 0.8 (flux); all three integrated power errors and the absolute balance must stay below 2e-8 W.

| Subdivisions per axis | Tetrahedra | Temperature RMS error (K) | Flux-vector RMS error (W/m²) | Temperature order | Flux order | Signed defect (W) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 48 | 0.766614 | 11.0627 | — | — | 2.86e-14 |
| 4 | 384 | 0.205439 | 5.91872 | 1.89979 | 0.90235 | 3.90e-13 |
| 8 | 3,072 | 0.0526435 | 3.02915 | 1.96438 | 0.96637 | 5.82e-12 |
| 16 | 24,576 | 0.0132662 | 1.52596 | 1.98850 | 0.98920 | 3.72e-11 |

## Serial process measurements

One run per grid on the same Linux workspace, AMD EPYC 9V74, without VTK export. Wall time includes load, mesh construction, assembly, solve, flux recovery and balance reassembly. RSS is the high-water mark returned by Linux wait4. Timings are observations of this host and run.

| Nodes | Tetrahedra | Time (s) | Peak RSS (KiB) | CG iterations | CG residual | Cold-face error (K) | Convective error (W) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 729 | 3,072 | 0.012967 | 10,496 | 71 | 3.16e-11 | 1.94e-11 | 1.95e-12 |
| 4,913 | 24,576 | 0.102274 | 19,100 | 137 | 3.17e-11 | 2.49e-11 | 3.01e-12 |
| 35,937 | 196,608 | 1.316938 | 137,280 | 253 | 3.84e-11 | 7.00e-11 | 6.08e-11 |

Reproduction commands and boundary conventions are in the [thermal project guide](../THERMAL_PROJECT.md). This evidence covers one stationary homogeneous box workflow. Experimental validation, general CAD, production scale and other physics require separate verification.
