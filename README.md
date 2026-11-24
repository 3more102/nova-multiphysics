# NOVA Multiphysics

**Engineering the Physics of Tomorrow.**

**NOVA Multiphysics** is an independent, C++20-based finite-element simulation project focused on building a verifiable numerical foundation for next-generation computer-aided engineering (CAE) software. Its long-term direction is a modular multiphysics platform spanning model definition, meshing, numerical solution, postprocessing, and a desktop engineering environment.

**Current release track:** `0.8.0` development SDK · **Platforms:** Windows and Linux · **Interface:** C++ library and command-line examples

> [!IMPORTANT]
> **Research and development software.** NOVA is an experimental numerical platform, **not** a certified, production-ready or general-purpose multiphysics product. Features have different levels of maturity, results require independent verification, and no claim of equivalence to commercial CAE packages is made.

[Get started](#quick-start) · [Capabilities](#implemented-capabilities) · [Examples](#reproducible-examples) · [Verification](#verification-and-quality) · [Architecture](#architecture-and-project-structure) · [Roadmap](docs/ROADMAP.md)

## Why NOVA?

NOVA is being developed around three engineering principles:

1. **Numerical credibility.** Analytical solutions, manufactured-solution studies, conservation diagnostics, and convergence checks take priority over feature counts.
2. **Transparent, reusable algorithms.** Finite elements, sparse linear algebra, time integration, and physics interfaces are implemented as separable C++ components with explicit assumptions.
3. **A foundation for growth.** The current FEM core and installable SDK are stepping stones toward broader physics, coupling strategies, model management, and visualization—not substitutes for those missing systems.

## Implemented capabilities

The following reflects functionality represented in the current source tree. **Experimental** indicates a limited-scope API or study, not a qualified engineering workflow.

| Area | Current implementation | Maturity |
| --- | --- | --- |
| **1D FEM** | Linear P1 Poisson; nonlinear P1 diffusion with Newton iterations and Armijo backtracking | Numerical prototypes |
| **2D scalar FEM** | Continuous triangular P1/P2 diffusion; spatially varying positive isotropic coefficients; Dirichlet, Neumann and Robin conditions | Numerical prototypes |
| **2D transient heat** | Consistent capacity matrix; Backward Euler and Crank–Nicolson; time-dependent loads and boundary data; adaptive timestep control and restart | Experimental |
| **3D steady diffusion / conduction** | Continuous P1 tetrahedra; isotropic or SPD tensor coefficients; material regions, heat sources, mixed boundaries and flux diagnostics | Experimental |
| **3D transient heat** | P1 tetrahedra; implicit time integration, adaptive steps, per-step energy accounting, time-varying boundary data, state checkpoints | Experimental |
| **3D electrothermal FEM** | Stationary DC conduction → Joule heating → thermal diffusion; temperature-dependent sigma(T)/k(T), segregated implicit transient coupling, adaptive step doubling and versioned step-boundary restart | Experimental C++ APIs |
| **Optics reference models** | Normal-incidence multilayer scattering, passive Drude dispersion, planar SPP and prescribed-intensity Kerr phase | Experimental analytical C++ API; separate from Maxwell FEM |
| **3D electromagnetics** | Lowest-order tetrahedral Nédélec edge elements, frequency-domain Maxwell system with complex material permittivity and PEC boundaries, field reconstruction and cell absorption; small lossless PEC cavity eigenmode reference solver | Experimental C++ APIs; no ports or PML |
| **3D electrostatics** | Scalar/tensor permittivity, electric potential, field/displacement and energy; electrode-charge and capacitance diagnostics; restricted floating-terminal Schur workflow | Experimental |
| **3D linear elasticity** | Small-strain isotropic P1 tetrahedral solver with fully prescribed nodal displacements, body forces, face tractions, and stress/strain recovery | Experimental C++ API |
| **Mesh infrastructure** | Validated 2D triangle/3D tetrahedron topology, quality metrics, 2D uniform refinement, 3D uniform and marked local refinement, field transfer | Developing |
| **3D error indicators / adaptation** | P1 cell residuals, interior flux jumps and selected natural-boundary residuals; restricted one-cycle solve–estimate–mark–refine–resolve workflow | Experimental |
| **Numerical solvers** | CSR sparse algebra; CG, BiCGSTAB, restarted GMRES; Jacobi and CG-only IC(0); Newton nonlinear solves | Core components |
| **Interoperability / output** | First-order Gmsh MSH 2.2/4.1 ASCII import; legacy VTK, ParaView PVD time series; installable CMake SDK | Supported development workflows |

**Important scope boundaries:** NOVA contains experimental dedicated electrothermal coupling, but not a unified general-purpose multiphysics coupling framework. The project does not currently provide an integrated desktop GUI, general CAD kernel, automatic production-grade 3D meshing, validated fluid or general-purpose EM solvers, or commercial engineering certification.

### Physics and numerical documentation

See [optics reference models and limitations](docs/OPTICS_FOUNDATION.md).

| Topic | Technical reference |
| --- | --- |
| Steady 3D conduction, materials, boundary groups and verification | [Steady conduction on tagged 3D solids](docs/STEADY_CONDUCTION_3D.md) |
| Transient 3D heat, time integrators and energy diagnostics | [Transient 3D heat](docs/TRANSIENT_HEAT_3D.md) |
| Imported-mesh transient CLI | [Gmsh transient heat](docs/GMSH_TRANSIENT_HEAT_3D.md) |
| 2D/3D temporal adaptation | [Adaptive heat 2D](docs/ADAPTIVE_HEAT_2D.md) · [Adaptive heat 3D](docs/ADAPTIVE_HEAT_3D.md) |
| Electrothermal adaptive stepping and restart | [Adaptive electrothermal integration](docs/ADAPTIVE_TRANSIENT_ELECTROTHERMAL_3D.md) · [Versioned checkpoint](docs/ADAPTIVE_ELECTROTHERMAL_CHECKPOINT_3D.md) |\n| Electrostatics, capacitance and floating conductors | [3D electrostatics](docs/ELECTROSTATICS_3D.md) |
| 3D mesh refinement and residual indicators | [Tetrahedral refinement](docs/UNIFORM_TETRA_REFINEMENT_3D.md) · [Residual error indicators](docs/RESIDUAL_ERROR_ESTIMATOR_3D.md) |
| Conventions and analytical verification | [Numerical conventions](docs/NUMERICAL_CONVENTIONS.md) · [3D spatial error norms](docs/VERIFICATION_3D_ERROR_NORMS.md) |

## Quick start

### Requirements

- A **C++20** compiler (GCC 11+, Clang 14+, or MSVC 2022 recommended).
- **CMake 3.21+** and **Ninja** for the included presets. Generator-independent builds require CMake 3.20+.
- **Python 3** is optional; some independent verification scripts are registered when Python is available.
- Supported development targets include **Windows and Linux**.

Clone the project, configure a Release build, compile, and run the registered tests:

~~~bash
git clone https://github.com/3more102/nova-multiphysics.git
cd nova-multiphysics

cmake --preset default
cmake --build --preset default
ctest --preset default
~~~

Alternatively, without CMake presets:

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
~~~

On Windows with a multi-configuration generator, executables are typically found in `build/Release/`; with the default Ninja preset, they are in `build/`. A successful compilation is not a substitute for running the tests.

## Reproducible examples

Run these commands from the repository root after building the examples.

### 1. Steady-state heat conduction through a layered 3D wall

The sample mesh has two material regions and tagged exterior surfaces. This run applies fixed temperature, convection, and two thermal conductivities, then exports a ParaView-readable VTK file.

~~~bash
./build/nova_gmsh_diffusion_3d \
  examples/gmsh_diffusion_3d/two_layer_wall.msh \
  --material 11 2 --material 22 6 \
  --temperature 101 400 \
  --convection 102 12 300 \
  --balance-report --vtk wall.vtk
~~~

For the **zero-volumetric-source** reference case, the analytical 1D solution gives a **340 K material-interface temperature** and **240 W total heat transfer** through the 1 m² wall. The reference is a verification case, not evidence of accuracy for arbitrary geometries or materials.

Inspect mesh quality without performing a solve:

~~~bash
./build/nova_gmsh_diffusion_3d \
  tests/data/unit_tetrahedron.msh \
  --check-mesh --min-mean-ratio 0.8
~~~

Read the [model conventions, SI units, quality checks and verification criteria](docs/STEADY_CONDUCTION_3D.md).

### 2. Transient heat in a 3D solid

~~~bash
./build/nova_gmsh_transient_heat_3d \
  examples/gmsh_diffusion_3d/two_layer_wall.msh \
  --material 11 2 --material 22 6 \
  --capacity 11 2000000 --capacity 22 3500000 \
  --initial 300 --dt 0.5 --steps 40 --scheme cn \
  --temperature 101 310 --convection 102 12 290 \
  --heat-source 22 120 \
  --vtk-series heated_wall.pvd \
  --energy-csv heated_wall_energy.csv
~~~

Open `heated_wall.pvd` in **ParaView** to visualize temperature over time. The CSV provides per-step discrete energy accounting. Consult the [transient CLI documentation](docs/GMSH_TRANSIENT_HEAT_3D.md) for exact boundary semantics, supported mesh formats and limitations.

### 3. Additional numerical demonstrations

~~~bash
./build/nova_poisson_1d
./build/nova_diffusion_2d
./build/nova_transient_heat_2d
./build/nova_adaptive_heat_2d
./build/nova_nonlinear_diffusion_1d
./build/nova_electrostatics_3d
~~~

The 3D elasticity and selected advanced electrostatic/adaptation workflows are currently exposed through **C++ APIs and test cases**, not dedicated end-user CLI applications.

## C++ SDK integration

NOVA's core can be installed and consumed using CMake package discovery:

~~~bash
cmake --install build --config Release --prefix ./nova-sdk
~~~

In a separate C++20 CMake project:

~~~cmake
find_package(NOVA 0.8 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE NOVA::Core)
~~~

Configure the consuming project with `CMAKE_PREFIX_PATH` pointing to the installation prefix. The package is **pre-1.0** and has **no guaranteed stable API or ABI**. For full external-consumer compilation and testing instructions, see [SDK installation and verification](docs/SDK_INSTALL.md).

## Verification and quality

NOVA's verification strategy combines multiple forms of evidence rather than treating a passing solver iteration as proof of physical correctness.

- **Analytical reference problems:** affine patches, layered conduction, capacitor fields and elastic displacement/stress cases.
- **Manufactured solutions:** controlled spatial and temporal convergence studies, including L2/H1 error measurements.
- **Conservation diagnostics:** discrete thermal energy budgets, weak-form reactions, flux jumps and selected face-integrated heat reconciliation.
- **Robustness tests:** invalid mesh topology, non-finite inputs, unsupported boundary configurations, incompatible Neumann data and solver breakdowns.
- **Regression automation:** CTest-registered C++ tests, CLI smoke tests and optional Python verification scripts.

Run the currently configured suite with:

~~~bash
ctest --preset default --output-on-failure
~~~

The command reports the results for **your build and environment**. No numerical test pass count, platform matrix or performance benchmark is asserted without a corresponding run or published artifact.

**Interpretation matters:** discrete conservation, algebraic convergence and raw P1 flux recovery do not establish continuum accuracy, mesh independence, experimental agreement or engineering certification. See the [verification plan](docs/VERIFICATION_PLAN.md).

## Architecture and project structure

NOVA separates numerical infrastructure from physics modules and integration boundaries:

~~~text
include/nova/      Public C++ API headers
src/linalg/        Sparse matrices, iterative linear solvers and preconditioners
src/fem/           Meshes, finite elements, assembly, time integration and indicators
src/physics/       Physics-facing interfaces and operators
src/io/            Mesh import, VTK/PVD output and checkpoints
examples/          Standalone numerical demonstrations and sample meshes
tests/             Unit, regression, CLI and package-consumer tests
verification/      Independent analytical and result-checking scripts
docs/              Formulations, assumptions, limitations and development roadmap
cmake/             Build configuration and SDK package infrastructure
~~~

**Present:** an independently built FEM-oriented C++ core, solver library, mesh I/O, physics prototypes and installable SDK.

**Planned—not yet an integrated product:** desktop UI, geometry/CAD services, generalized material and field registries, segregated/monolithic multiphysics coupling, broader physics, advanced postprocessing and automation.

Explore the [architecture design](docs/ARCHITECTURE.md) and [milestone roadmap](docs/ROADMAP.md). Milestones require documented acceptance tests and verification evidence; an implemented prototype alone does not signify product readiness.

## Engineering policies

- Document **units, signs, boundary conditions, assumptions and unsupported cases**.
- Report solver convergence and failures explicitly; do not silently substitute algorithms.
- Keep independently reproducible numerical benchmarks and regression tests.
- Review third-party dependencies for security, maintenance and future commercial licensing.
- Separate **verified implementation**, **experimental research** and **planned product features** in all external claims.

## Licensing and usage

**A public repository does not automatically grant permission to reuse or redistribute its contents.** NOVA's project-wide license and commercial distribution terms have not been finalized. Do not assume open-source or commercial redistribution rights, or redistribute the prototype commercially, without explicit authorization and a completed dependency/license review.

## Disclaimer

NOVA Multiphysics is ongoing experimental software research. It is **not safety-certified, production-qualified or validated for consequential engineering design decisions**. Independently check modeling assumptions, boundary conditions, spatial and temporal convergence, solver diagnostics and applicable standards before using any computed result.

---

**NOVA Multiphysics** · *Independent numerical engineering. Verification first.*
