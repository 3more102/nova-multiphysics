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

The verification code integrates the exported P1 field over every tetrahedron using its four nodal values and the tetrahedron's geometric volume. It checks that the summed mesh volume equals the analytic cube volume within 2e-12 m³; the mean calculation does not rely on equal element volumes. Independent acceptance criteria are a finest mean error below 0.025 K, lower bounds on both refinement orders above 1.65, all prescribed boundary nodes within 2e-8 K of 300 K, and each integrated power and balance within 2e-8 W of its reference. The steady generated and removed powers are both $QL^3=30$ W.

| Mesh | FEM mean rise (K) | Upper bound on error (K) | Lower bound on order | Signed heat defect (W) |
| --- | ---: | ---: | ---: | ---: |
| 4 × 4 × 4 | 0.426815257353 | 0.17824038 | — | 1.66e-12 |
| 8 × 8 × 8 | 0.552558507149 | 0.05249713 | 1.76350 | 9.03e-12 |
| 16 × 16 × 16 | 0.591197174132 | 0.01385846 | 1.92143 | 4.37e-11 |

The checker also verifies the full nodal field's inversion symmetry through the cube center (maximum reflected-node discrepancy below 2e-8 K), confirms that the unique center node contains the maximum nodal temperature, and checks that the runner's reported maximum matches this center probe. These are symmetry and peak-location checks, not an absolute continuum bound on the peak value.

## SciPy implementation and independent discretization checks

The first SciPy reference is a node-centred seven-point finite-difference solve on the same Cartesian vertices. A type-I discrete sine transform diagonalizes its operator. For this particular uniform cube, homogeneous conductivity and uniform source, the six-tetrahedron P1 Galerkin assembly reduces algebraically to the same seven-point nodal stencil. The SciPy code uses a separate implementation, matrix construction and solver, so the full-field comparison checks NOVA's assembly, solve and VTK output; it is not an independent spatial discretization. At 16 subdivisions, all 4,913 NOVA nodal temperatures are compared at identical coordinates, with maximum and RMS differences each below 2e-8 K.

| Subdivisions per axis | Nodal SciPy mean rise (K) | Upper bound on continuum mean error (K) | Center rise (K) |
| ---: | ---: | ---: | ---: |
| 16 | 0.591197174133 | 0.01385846 | 1.676429964553 |
| 32 | 0.601530120040 | 0.00352551 | 1.683880381680 |
| 64 | 0.604169110303 | 0.00088652 | 1.685757768523 |

The observed order from successive center-value differences between 16, 32 and 64 subdivisions is 1.9886. At 16 subdivisions, the nodal SciPy mean and center values agree with NOVA within the 0.025 K comparison allowance. The continuum Fourier-series bracket independently bounds the volume-mean error; the center-value sequence is a refinement check, not a rigorous continuum error bound for the peak temperature.

The second SciPy reference is a cell-centred finite-volume discretization on equal-volume control volumes. Dirichlet faces are half a cell from the first and last unknowns, giving a one-dimensional operator with diagonal 3 at the two boundary cells, diagonal 2 in the interior, and -1 on adjacent cells. The orthonormal type-II discrete sine transform diagonalizes this operator; the three-dimensional eigenvalues are sums of the one-dimensional eigenvalues. This is a genuinely different spatial discretization from NOVA's vertex-based P1 tetrahedra. Its volume mean is the arithmetic mean of the cell-centre values. [SciPy multidimensional sine transform](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.dstn.html).

| Subdivisions per axis | SciPy cell-FV mean rise (K) | Upper bound on continuum error (K) | Lower bound on order |
| ---: | ---: | ---: | ---: |
| 16 | 0.619850563835 | 0.014796136657 | — |
| 32 | 0.608804622953 | 0.003750195774 | 1.98007 |
| 64 | 0.605996464954 | 0.000942037775 | 1.99265 |

The cell-centred FV mean approaches the continuum bracket from above, while the NOVA P1 and nodal SciPy means approach it from below. The checker requires successive lower bounds on the finite-volume refinement order to exceed 1.8 and the finest upper error bound to be below 0.002 K.

Recorded using Debian 13 x86-64, GCC 14.2 Release, Python 3.12.14 and SciPy 1.17.0 with solver commit `3ab9afa446ea92392e8910fd42b90ce24d221eea`. The continuum comparison uses Python's standard library and runs in CTest. Ubuntu CI also runs both SciPy checks against the installed runner and installed project. To repeat the external calculations:

~~~sh
python verification/verify_thermal_project.py cube build/nova-run examples/projects/internally_heated_cube.nova --scipy
~~~

The evidence covers a homogeneous stationary cube, its global mean, nodal implementation agreement, inversion symmetry, center-peak location and integrated boundary reactions. It does not establish experimental accuracy or local peak-temperature tolerances for other geometries.
