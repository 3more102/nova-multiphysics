# Self-consistent 3D thermo-optic feedback

NOVA's one-way Maxwell-to-thermal solver transfers the exactly integrated dielectric loss from each Nédélec tetrahedron into the corresponding P1 thermal cell. This extension closes the stationary feedback loop by updating the optical permittivity from the cell-average temperature and repeating both solves.

For cell c, the supported constitutive law is

$$
\varepsilon_{r,c}(T_c)=\varepsilon_{r,c}^{(0)}
+\alpha_c(T_c-T_{\mathrm{ref}}),
$$

where $\alpha_c=d\varepsilon_{r,c}/dT$ is a complex, constant per-kelvin coefficient and $T_c$ is the arithmetic mean of the four P1 vertex temperatures. Both the real thermo-optic shift and the temperature dependence of dielectric loss are supported. Updated materials must remain finite and passive for the exp(-i omega t) convention.

## Iteration

1. Solve the base stationary thermal problem without optical heating to obtain an initial temperature field.
2. Evaluate cell permittivities from the current P1 temperature field.
3. Assemble and solve the complex 3D Maxwell system.
4. Integrate dielectric absorption per tetrahedron and transfer the cell power conservatively into the thermal source.
5. Solve the stationary P1 heat problem and under-relax the temperature iterate.
6. Stop only when both the maximum nodal temperature residual and relative change in total absorbed optical power satisfy their tolerances.

The result returns the final one-way coupled solve, the optical materials used for its Maxwell solve, and per-iteration residual/power history. Reaching the iteration limit returns the last finite result with `converged=false`; invalid inputs or failed inner linear solves throw.

## Scope and limitations

This is a stationary, one-way frequency-domain Maxwell / steady thermal fixed-point model with cellwise linear temperature-dependent isotropic permittivity. It preserves the existing thermal material and source models and does not mutate the input problem. It is not a transient thermo-optic solver, a temperature-dependent permeability/conductivity model, a Kerr/nonlinear-optical solver, or a port/PML implementation. No experimental validation or general convergence guarantee is implied; choose relaxation and tolerances for the specific material and geometry.


## Manufactured nonzero-temperature reference

A second regression isolates the temperature-dependent material update from the free-temperature fixed-point test. On a unit cube, all P1 thermal vertices are prescribed at 320 K, with $T_{\mathrm{ref}}=300$ K, $\varepsilon_r^{(0)}=2+0.6i$, and $d\varepsilon_r/dT=10^{-3}+2\times10^{-4}i$ K$^{-1}$. The expected uniform optical permittivity is therefore $2.02+0.604i$.

The cell source is manufactured as $F=-k_0^2\varepsilon_r(320)E_0$ for $E_0=(100,0,0)$, making the exact Maxwell field constant. The independent absorbed-power density is $q=\tfrac12 k_0c\varepsilon_0\operatorname{Im}(\varepsilon_r)|E_0|^2=18.4375924846$ W/m$^3$; the unit cube therefore absorbs 18.4375924846 W. The regression checks the updated permittivity in every cell, each integrated cell power, the global thermal reaction and zero energy defect. This complements the separate free-temperature nonzero-slope convergence test.

## Coupled equilibrium with four free thermal nodes

A regular tetrahedron with vertices (1,1,1), (1,−1,−1), (−1,1,−1) and
(−1,−1,1) m has volume V = 8/3 m³ and total area A = 8√3 m². Each face
exchanges heat with a 300 K ambient at h = 0.2 W/(m²·K). There is no
prescribed nodal temperature or background heating. Symmetry gives each P1
vertex one quarter of both the optical heat load and the convective heat
removal. At a uniform nodal temperature the conduction term vanishes.

For the natural Maxwell boundary and a constant x-directed source
−k₀²(2 + 0.6i) × 100, the electric field has zero curl and its only nonzero
component is Eₓ(T) = 100(2 + 0.6i)/ε(T). Here k₀ = 2.3 m⁻¹ and
ε(T) = 2 + 0.6i + (0.02 + 0.003i)(T − 300 K), with the slope per kelvin.

Solve the independent scalar equation
**½ k₀ c ε₀ Im(ε(T)) |Eₓ(T)|² V = h A (T − 300 K)**. This calculation
uses neither the NOVA assembly nor its iteration. Reproduce the reference
with SciPy:

~~~python
from math import sqrt
from scipy.optimize import brentq

epsilon_ref = 2 + 0.6j
def field(T):
    epsilon = epsilon_ref + (0.02 + 0.003j) * (T - 300)
    return epsilon, 100 * epsilon_ref / epsilon

def power(T):
    epsilon, ex = field(T)
    return 0.5 * 2.3 * 299792458 * 8.8541878128e-12 * epsilon.imag * abs(ex)**2 * (8/3)

T = brentq(lambda T: power(T) - 0.2 * 8 * sqrt(3) * (T - 300), 300, 340)
print(T, power(T), field(T)[1])
~~~

| Quantity | Reference |
| --- | ---: |
| All four nodal temperatures | 314.5552553103854 K |
| Absorbed power and convective removal | 40.3367067435588 W |
| Complex x-directed electric field | 87.72733374282475 + 1.54201504058201i V/m |

The regression checks the equilibrium with relaxation factors 0.5 and 1,
including all six complex edge integrals and thermal power balance. Temperature
tolerance against the reference is 2e-8 K, integrated-power tolerance is
2e-7 W, and complex edge-integral error tolerance is 2e-7 V. The uniform
temperature is an exact result for this one-element P1 model. Continuum
temperature accuracy and spatial convergence require a separate benchmark.
