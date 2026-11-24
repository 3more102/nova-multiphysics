# Heat-induced displacement and stress on a 3D Gmsh solid

Both Gmsh heat executables can solve quasistatic mechanical equilibrium from the temperature they compute. The transient program can evaluate mechanics on the Dirichlet-projected initial state and every completed Backward Euler or Crank–Nicolson step. Temperature and displacement share the original P1 mesh vertex numbering. Region coefficients, recovered flux, strain and mean stress follow the imported tetrahedron ordering.

## Run the examples

From a Release build, reproduce the stationary block, mixed-boundary layered wall, rotated-tensor diffusion, free expansion and transient cooling studies in one command. The verifier retains readable VTK frames, PVD collections and thermal energy CSV files when given an empty output directory.

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
python verification/verify_thermoelastic_3d.py build/nova_gmsh_diffusion_3d build/nova_gmsh_transient_heat_3d --output results
~~~

For a multi-configuration Windows build, pass the two executables in `build/Release/`, with their `.exe` extensions. To run a single cooling study from the repository root, use the following command.

~~~sh
./build/nova_gmsh_transient_heat_3d examples/thermoelastic_3d/tetrahedron.msh \
  --material solid 1 --capacity solid 1 --initial 310 --temperature base 300 \
  --dt 0.025 --steps 8 --scheme cn --elastic solid 100 0 0.01 300 \
  --clamp base --vtk-series cooling.pvd --energy-csv cooling.csv
~~~

Open `cooling.pvd` in ParaView. Color by point temperature or cell mean stress. Warp by the point displacement vector to view deformation. The geometry stored on disk uses the original coordinates.

| Argument | Meaning |
| --- | --- |
| `--elastic VOL E nu alpha T_ref` | Supply one law for every physical volume. Young's modulus is in Pa, thermal expansion in 1/K and the reference temperature uses the same scale as the computed field. Each region can have a different reference. |
| `--clamp SURF` | Prescribe all three displacement components to zero on a physical surface. |
| `--fix-component SURF x\|y\|z value` | Prescribe a global Cartesian component in metres. Compatible values on intersecting surfaces share a DOF. |
| `--elastic-solver jacobi\|ic0 rtol iterations` | Configure mechanical CG. The CLI defaults to Jacobi, relative tolerance 1e-11, absolute tolerance 1e-14 and 1,000 iterations. Failure aborts publication of the affected frame. |
| `--vtk` or no visualization option | Stationary and final-only transient runs compute a single mechanical equilibrium. The CLI reports residual, iteration count and stored elastic energy. |
| `--vtk-series` | Transient runs compute initial and accepted-frame equilibria. Mesh geometry, stiffness, prescribed-column corrections and eligible IC(0) factors are reused. |

A steady coupled study needs a thermal Dirichlet datum or positive Robin transfer on every component. A free additive temperature constant would change the predicted expansion. A transient study obtains its temperature datum from the initial condition. Each face-connected elastic solid needs supports of rank six on its rigid translations and rotations. The rank test uses centered, scaled anchor coordinates and a 1e-12 numerical rank threshold. Bodies connected only at a vertex or edge must each satisfy that support check.

## Discretization and output

For isotropic small-strain elasticity, let `s = alpha*(T - T_ref)`, `eps = sym(grad u)`, and `B = E/[3(1-2 nu)]` be the scalar thermal strain, total strain, and bulk modulus. The constitutive equation is `sigma = C:(eps - s I)`. Known thermal strain contributes a positive load `integral 3 B s grad(N_a) dV` to the displacement equation. Quasistatic mechanics receives temperature after the heat solve. The thermal energy budget accounts for the thermal equation alone. [Continuum formulation, J. Bleyer](https://bleyerj.github.io/comet-fenicsx/tours/multiphysics/thermoelasticity_weak/thermoelasticity_weak.html)

On a straight P1 tetrahedron, displacement strain and basis gradients are constant. With constant material properties per cell and four nodal thermal strains `s_i`, the load and cell-mean stress use `m = sum(s_i)/4`. Exact discrete elastic energy also needs the variance across the cell.

~~~text
v = sum_i (s_i - m)^2 / 20
mean_stress = C:(eps - m I)
U_cell = volume/2 * [(eps - m I):mean_stress + 9 B v]
~~~

The tetrahedral consistent mass moments give this variance directly. Omitting it leaves the displacement unchanged and underestimates the energy whenever the nodal thermal strains differ. In the single-tetrahedron test with tip temperature 310 K, the centroid-only energy is 0.01041666667 J and the integrated energy is 0.01979166667 J.

The public `PreparedElasticity3D` owns the geometry and constrained operator. It accepts P1 nodal temperature and region laws through `solve_temperature`, or cellwise isotropic eigenstrain through `solve_cell_eigenstrain`. The original one-shot elasticity entry point uses the same assembly. Sparse assembly gathers one row of 3×3 nodal blocks at a time, eliminates constrained columns, then transfers sorted rows directly into CSR. Warm starts retain only a converged iterate. A repeated solve uses the same free-force residual scale as its corresponding cold start.

Coupled VTK files contain point temperature and displacement, cell material ID, heat flux, temperature gradient, total strain and mean stress. Available generation data use a cell scalar. The stress and strain arrays are symmetric 3×3 tensors with tensor shear. Global field data store the mechanical CG iterations, residual norm and integrated elastic energy. [VTK legacy format specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtk_legacy_file_format.html)

## Numerical evidence

The measured cases below use the Release build on Linux x86-64. All exported tensor components and field associations are checked by the end-to-end verifier.

| Model | Analytical reference | Measured |
| --- | --- | --- |
| Two layers with both end temperatures fixed | Interface 325 K, axial flux 300 W/m² | 325 K, 300 W/m², elastic energy 0.398125 J |
| Layered wall with Dirichlet, Neumann and Robin faces | Interface 340 K, axial flux 240 W/m² | 340 K, 240 W/m², elastic energy 0.910000 J |
| Rotated SPD conductivity, three prescribed outward fluxes | T = 300 + 10z K, q = (−5, 2.5, −40) W/m² | Tip 309.99999999999994 K, mechanical energy 0.01979166666666644 J |
| Uniform heating on three Cartesian symmetry planes | u = 0.05 (x,y,z), zero stress | Maximum absolute stress 1.39e-15 Pa |
| Cooling with the base held at 300 K | Consistent P1 mass and stiffness give a' = −10a for the free-node excess | At 0.2 s and dt = 0.025 s the BE tip is 301.677721600 K and CN is 301.339196312 K. Maximum step heat imbalance is 5.33e-15 J. |

Absolute tolerances for the small fixtures are 3e-10 K for temperatures, 2e-8 W/m² for wall flux, 2e-9 for single-tetrahedron flux/stress components, and 3e-12 J for single-tetrahedron energy. These bounds allow accumulation and iterative solve error above double roundoff while resolving the sign, tensor ordering and variance contribution. Cooling is compared to the exact semidiscrete recurrence at each exported step. Temporal truncation is measured separately against `300 + 10 exp(−10t)` at the tip.

| dt (s) | BE tip error (K) | CN tip error (K) |
| ---: | ---: | ---: |
| 0.025 | 0.324368768 | 0.014156521 |
| 0.0125 | 0.165653698 | 0.003528033 |
| 0.00625 | 0.083704056 | 0.000881319 |

The observed orders on the last halving are 0.985 and 2.001, respectively. The automated lower bounds are 0.90 and 1.95. A separate manufactured cube has `T=300+x`, `alpha=1` and stress-free displacement `((x²-y²-z²)/2, xy, xz)`. With exact boundary displacements, the P1 mechanical solution has free interior DOFs. NOVA integrates displacement errors against the polynomial with its independent Duffy quadrature implementation.

| Subdivisions per axis | Displacement L2 error | Displacement H1 seminorm error | Stored energy (J), exact 0 |
| ---: | ---: | ---: | ---: |
| 2 | 0.04564354646 | 0.4787135539 | 0.8750000 |
| 4 | 0.01141088661 | 0.2393567769 | 0.2187500 |
| 8 | 0.002852721655 | 0.1196783885 | 0.0546875 |

## Performance

The comparison uses upstream revision `16ef7f3dcd991e37623948c115856954a00dbd20` and the same benchmark driver, Linux/GCC 14.2.0 Release (-O3), AMD EPYC 9V74, single-threaded solves, and three fresh processes per case. Original and prepared runs alternate order. Each imports the same two-region Gmsh mesh, solves five CN steps, evaluates six mechanical equilibria with Jacobi CG, and writes one final thermal VTK. The spatially uniform cooling field makes the old cellwise temperature input and the new P1 transfer physically equivalent. The driver checks temperature against the exact ramp and compares the two executables' displacement norms, stress norms and summed energies at relative tolerance 2e-8.

| Tetrahedra | Thermal vertices | Original time (s) | Prepared time (s) | Original / prepared peak RSS (MiB) | Full CLI, six coupled VTK frames (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 384 | 125 | 0.0139 | 0.00565 | 4.56 / 2.68 | 0.0268 |
| 3,072 | 729 | 0.1385 | 0.0660 | 20.68 / 7.16 | 0.1841 |
| 24,576 | 4,913 | 1.9300 | 0.8204 | 149.10 / 41.41 | 1.6403 |
| 82,944 | 15,625 | 8.1196 | 4.0002 | 491.70 / 132.84 | 8.2447 |

Times and RSS are medians. The independent full CLI column includes writing all thermal and mechanical arrays for each frame. Its peak RSS at the largest size was 132.23 MiB. The original SDK had no such multi-field CLI, so its time comparison uses the common one-file workload.

On 82,944 tetrahedra, the prepared driver measured 0.164 s to import, 0.0448 s for mechanical topology, 0.0060 s for element/constraint/static-load setup and 0.0480 s for fused element-block and global CSR stiffness assembly. Across six mechanical fields, loads took 0.0263 s, CG took 3.0963 s, and recovery took 0.0230 s. Thermal integration excluding the mechanical callbacks took 0.2807 s, and final flux recovery plus thermal VTK writing took 0.1729 s. Phases are independently medianed. Peak RSS comes from the kernel's child-process high-water mark.

Reproduce the prepared workload and CLI measurements with the following commands. GNU time is optional. For a historical comparison, compile the same driver against the baseline checkout and pass its executable with `--reference`.

~~~sh
cmake --build build --config Release --target nova_benchmark_thermoelastic_3d
git worktree add --detach ../nova-16ef7f3 16ef7f3
cmake -S ../nova-16ef7f3 -B ../nova-baseline-build -DCMAKE_BUILD_TYPE=Release
cmake --build ../nova-baseline-build --target nova_core --parallel
c++ -std=c++20 -O3 -DNDEBUG -DNOVA_BENCH_REFERENCE -I../nova-16ef7f3/include \
  verification/benchmark_thermoelastic_3d.cpp ../nova-baseline-build/libnova_core.a -o ../nova-reference
python verification/benchmark_thermoelastic_3d.py build/nova_benchmark_thermoelastic_3d \
  build/nova_gmsh_transient_heat_3d --reference ../nova-reference --time-exe /usr/bin/time \
  --repeat 3 --sizes 4 8 16 24
~~~

## Scope

The model assumes infinitesimal strain, linear isotropic elasticity, cellwise constant expansion coefficients, a fixed conforming tetrahedral mesh and quasistatic response. Conductivity may be anisotropic. Stress is reported as a cell mean, and the consistent energy includes within-cell thermal variation. Deformation-to-heat feedback, inertial mechanical response, contact, plasticity, near-incompressible mixed elements and adaptive coupled remeshing are outside this implementation.
