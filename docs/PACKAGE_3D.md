# Coupled transient package analysis

The installed C++ interface is `nova/physics/package_3d.hpp`. A prepared model owns one straight P1 tetrahedral mesh, a region-law table, nodal terminal and temperature constraints, heat exchange on exterior triangles, and mechanical supports/loads. Each conducting component must have a prescribed voltage. Zero-conductivity regions generate zero current and Joule heat. At electrically inactive vertices the stored zero potential is a display convention identified by the active mask.

The transport unknown is ordered as all nodal temperatures followed by all nodal potentials. Current continuity and heat balance share one sparse residual. At cell-average temperature each scalar property is its reference value times `exp(slope*(T-Tref))`. Temperature changes electrical conductivity, thermal conductivity and Young's modulus independently in each region. Capacity is integrated into a material enthalpy. The consistent tetrahedral mass matrix acts on the P1 interpolant of the nodal enthalpies. This preserves a discrete first law for nonlinear capacity. Surface convection uses a degree-two rule. Diffuse-gray radiation to a prescribed large isothermal enclosure uses a face rule exact through degree five. Perfect thermal and mechanical bonding is assumed.

Damped Newton uses all four analytic transport tangent blocks with GMRES or BiCGSTAB. A segregated strategy alternates the voltage solve and a damped temperature solve and checks the complete residual. The scales for volts and kelvin are configurable. Quasistatic small-strain mechanics is the downstream block of a triangular system. The stiffness is reused for constant modulus and rebuilt at the accepted temperature for a variable modulus. Displacements do not change geometry or Joule heating. Mechanical work and thermoelastic damping are outside the heat budget.

Backward Euler, Crank–Nicolson and variable-step BDF2 share the residual. For timestep h, preceding step k and r=h/k, BDF2 uses (1+2r)/[(1+r)h] times the newest enthalpy increment minus r/[(1+r)k] times the preceding increment. Startup and the first step after a configured event use Crank–Nicolson for BDF2. Growth is capped at two. Endpoint averaging in CN applies to spatial heat loads, including Joule power. Electrical equilibrium always uses the current endpoint.

Adaptive BE accepts two half steps and compares them with one full step. Adaptive CN/BDF2 compare their endpoint temperature with an embedded BE solution. This difference is a first-order-pair indicator. It is not an exact high-order local error estimate. Only a complete accepted transport and mechanical equilibrium updates user state. The diagnostic heat residual is the signed sum of unconstrained physical thermal rows in watts. For BDF2 this uses the discrete enthalpy derivative and should not be interpreted as a simple one-step energy increment. Voltage reactions give terminal currents and terminal electrical power. Sources and convective/radiative power are reported separately.

Variable-step formulation reference. [Li and Liao, SIAM Journal on Numerical Analysis](https://doi.org/10.1137/21M1462398).

## Run a package analysis

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
./build/nova_package_3d examples/package_3d/power_module.nova package-results
python verification/verify_package_3d.py build/nova_package_3d --independent
~~~

The final command additionally requires NumPy and SciPy. The solver and default CTest verification have no such dependency. On a multi-config generator use the executable under `build/Release`. Open `package-results/package.pvd` in ParaView. Individual VTK files contain nodal temperature, potential, electrical activity, displacement and support reaction, and cell material tag, current, heat flux, Joule density, conductivity scale, Young's modulus, strain, stress, principal stresses and von Mises stress. CSV histories contain accepted-step diagnostics, volume-averaged and extremal temperatures, interpolated probes and per-terminal currents. The JSON summary records mesh size, rejection count and residual extrema.

### Model version 1

The first two records must be `NOVA_PACKAGE 1` and `mesh path`. Paths resolve relative to the model. Comments begin with a hash. All quantities are SI. Gmsh physical group names without spaces and numeric tags are supported. The tetrahedral mesh can be ASCII Gmsh 2.2 or 4.1. Each physical volume requires one material.

| Record | Arguments, in order |
| --- | --- |
| material | volume, reference K, electrical conductivity and log slope, thermal conductivity and log slope, volumetric heat capacity and log slope, Young's modulus and log slope, Poisson ratio, isotropic expansion per kelvin |
| tensor | electrical or thermal, previously defined volume, dimensionless SPD components xx yy zz xy xz yz. The scalar reference law multiplies this tensor. |
| curve | identifier, one or more time/value pairs. Linear interpolation with constant endpoint extrapolation. Every knot in a varying curve is a time event. |
| voltage / temperature | exterior physical group, previously defined curve |
| source | volume, additive heat generation curve in W/m³ |
| cool | exterior group, convection coefficient, emissivity, ambient K curve, optional outward flux curve |
| clamp | exterior group with zero displacement |
| initial | uniform initial temperature in kelvin |
| study | be, cn or bdf2, end time, initial step, minimum step, maximum step |
| solver | monolithic or segregated, scaled absolute residual tolerance, relative residual tolerance, iteration limit |
| linear | gmres or bicgstab, none/jacobi/ilu0, maximum linear iterations |
| adapt | absolute Kelvin and relative temperature tolerances. Both zero disables error control. |
| scales / relaxation | characteristic kelvin and volts / segregated relaxation factor |
| probe | identifier and x y z. Interpolation uses the containing P1 tetrahedron. Voltage in an insulating cell is left blank. |
| frames | write every Nth accepted VTK frame and the final frame. Zero disables VTK. CSV and checkpoints continue. |

An accepted event or final time writes a new checkpoint. Its record contains the exact mesh, two accepted temperature levels, potential and proposed step. A fingerprint covers the original model and mesh bytes including units encoded by the model, material and boundary tags, and solver settings. Editing either file changes its identity. This detects accidental mismatch and provides no cryptographic authentication. Restart into a fresh output directory with `--restart path/to/step_NNNNNN.chk`. Use `--until seconds` for an earlier stopping point. Stopping at an existing accepted event preserves the original step sequence. An inserted time boundary may change that sequence.

### Power-module geometry and assumptions

The supplied 4 mm × 2 mm coupon has a 0.4 mm isolated metal heat spreader, a 0.1 mm effective thermal bond, a 0.3 mm ceramic substrate, a 0.2 mm conducting semiconductor die and 0.1 mm top contacts separated by passivation. The voltage spans two top contacts and varies over 80 ms. The spreader has a zero reference and no conductive path to the die. The base is clamped and cooled. The exposed passivation exchanges heat with ambient by convection and radiation. Other exterior surfaces are insulated. The finite bond thickness represents a distributed through-thickness resistance; interfaces themselves are perfectly bonded. All tabulated material data, slopes, voltages and boundary coefficients are illustrative engineering inputs. No device measurement, doping calibration, junction physics, interfacial slip or failure criterion is assumed.

Regenerate the mesh with `python examples/package_3d/generate_mesh.py case.msh --level 2`. The reference solver independently builds conduction, elastic and enthalpy-mass matrices in SciPy and computes a time-continuous thermal reference on the same spatial discretization. Independent geometry, experimental material identification and component testing would be needed for product validation.
