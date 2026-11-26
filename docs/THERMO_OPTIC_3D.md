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


## Spatially varying manufactured reference

The affine-temperature regression prescribes all thermal vertices on the unit cube with $T(x,y,z)=300+10x+20y+30z$ K (coordinates in metres). For a P1 tetrahedron, its exact volume-average temperature is the affine field evaluated at the geometric centroid. With the same reference permittivity, slope and optical amplitude as the uniform reference above, the manufactured cell source is $F_c=-k_0^2\varepsilon_{r,c}E_0$. This supports the exact spatially constant optical field even though cell permittivities differ.

The local target is $q_c=\tfrac12 k_0c\varepsilon_0\operatorname{Im}(\varepsilon_{r,c})|E_0|^2$. Because the cube-average temperature is 330 K and the constitutive law is affine, the total absorbed power is **18.4986441154 W**. The test checks individual permittivities, integrated cell absorption, a nonuniform density, prescribed nodal temperatures and the global reaction balance. It detects unintended uniform material averaging or incorrect cell-to-cell loss transfer. This is a manufactured discrete consistency reference, not a free-temperature or experimentally validated optical model.

A separate negative regression starts with passive reference permittivity $2+0.6i$ and an imaginary slope of $-0.1$ K$^{-1}$ at 320 K. The updated imaginary part is $-1.4$; the passive-only solver must reject this update with an invalid-argument exception before the optical solve.
