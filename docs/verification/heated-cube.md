# Uniformly heated 3D cube

The native project `internally_heated_cube.nova` uses a cube of side length 0.5 m, isotropic conductivity 2 W/(m·K), volumetric generation 240 W/m³ and a fixed temperature of 300 K on all six faces. All quantities are synthetic verification inputs. Interior temperature varies along all three coordinates. The physical observables are volume-mean temperature rise, fixed-boundary heat removal and total generated power.

## Continuum reference

Define the rise above the boundary temperature as $u$. The boundary value problem is $-k\nabla^2u=Q$ with $u=0$ on all six faces. Expanding the constant source in the three Dirichlet sine bases and integrating the solution over the cube gives

$$
\bar u=\frac{512QL^2}{k\pi^8}
\sum_{a,b,c\;\mathrm{odd}>0}
\frac{1}{a^2b^2c^2(a^2+b^2+c^2)}.
$$

For a cutoff $N$ on each index every term is positive. Bound the omitted terms with $a>N$ by replacing $a^2+b^2+c^2$ with $a^2$. Use $\sum_{a>N}a^{-4}\le1/(3N^3)$, $\sum_{j\;\mathrm{odd}>0}j^{-2}=\pi^2/8$, and sum the three axis tail bounds. The truncation error is at most $8QL^2/(k\pi^4N^3)$. At $N=127$ the resulting mean lies between 0.605054427178 and 0.605055629998 K (displayed endpoints rounded). The bound concerns truncation of the continuum series. [Dirichlet eigenfunction expansion for the Poisson equation, UCLA](https://www.math.ucla.edu/~yanovsky/handbooks/PDEs.pdf).

The verification code integrates the exported P1 field over every tetrahedron using its four nodal values. All tetrahedra in this structured cube have equal volume. Independent acceptance criteria are a finest mean error below 0.025 K, lower bounds on both refinement orders above 1.65, all prescribed boundary nodes within 2e-8 K of 300 K, and each integrated power and balance within 2e-8 W of its reference. The steady generated and removed powers are both $QL^3=30$ W.

| Mesh | FEM mean rise (K) | Upper bound on error (K) | Lower bound on order | Signed heat defect (W) |
| --- | ---: | ---: | ---: | ---: |
| 4 × 4 × 4 | 0.426815257353 | 0.17824038 | — | 1.66e-12 |
| 8 × 8 × 8 | 0.552558507149 | 0.05249713 | 1.76350 | 9.03e-12 |
| 16 × 16 × 16 | 0.591197174132 | 0.01385846 | 1.92143 | 4.37e-11 |

The checker also verifies the full nodal field's inversion symmetry through the cube center (maximum reflected-node discrepancy below 2e-8 K), confirms that the unique center node contains the maximum nodal temperature, and checks that the runner's reported maximum matches this center probe. These are symmetry and peak-location checks, not an absolute continuum bound on the peak value.

## Independent finite-difference calculation

The optional SciPy checker solves the seven-point finite-difference equations on a separate Cartesian grid. A type-I discrete sine transform diagonalizes the operator and an inverse transform returns its temperature rise. The grid mean uses tensor-product trapezoidal integration, with zero rise on the boundary. It shares no NOVA mesh, matrix, linear solver or output fields. The finest error allowance is 0.002 K and successive orders must exceed 1.8. [SciPy multidimensional sine transform](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dstn.html).

| Subdivisions per axis | SciPy mean rise (K) | Upper bound on continuum error (K) |
| ---: | ---: | ---: |
| 16 | 0.591197174133 | 0.01385846 |
| 32 | 0.601530120040 | 0.00352551 |
| 64 | 0.604169110303 | 0.00088652 |

The same independent finite-difference solution also gives the center temperature rise:

| Subdivisions per axis | Center rise (K) |
| ---: | ---: |
| 16 | 1.676429964553 |
| 32 | 1.683880381680 |
| 64 | 1.685757768523 |

The observed order from successive center-value differences between 16, 32 and 64 subdivisions is 1.9886. At 16 subdivisions, the independent center value agrees with NOVA's center nodal result within the 0.025 K comparison allowance. This supports the peak-location and local refinement checks, but is not a rigorous continuum error bound for the peak temperature.

Recorded using Debian 13 x86-64, GCC 14.2 Release, Python 3.12.14 and SciPy 1.17.0 with solver commit `3ab9afa446ea92392e8910fd42b90ce24d221eea`. The continuum comparison uses Python's standard library and runs in CTest. Ubuntu CI also runs the SciPy comparison against the installed runner and installed project. To repeat the external calculation:

~~~sh
python verification/verify_thermal_project.py cube build/nova-run examples/projects/internally_heated_cube.nova --scipy
~~~

The evidence covers a homogeneous stationary cube, its global mean, nodal inversion symmetry, center-peak location, and integrated boundary reactions. It does not establish experimental accuracy or local peak-temperature tolerances for other geometries.
