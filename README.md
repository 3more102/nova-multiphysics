# NOVA Multiphysics

**Engineering the Physics of Tomorrow.**

NOVA is an independent, C++20 multiphysics simulation platform under development. Its long-term objective is to provide an integrated environment for geometry, meshing, finite-element modeling, numerical solution, multiphysics coupling, and scientific visualization.

> **Development status:** Early engineering prototype. Not production-ready and not yet validated for engineering design decisions.

## Product direction

NOVA is being developed as a commercial, general-purpose computer-aided engineering (CAE) platform with an independently developed numerical core.

## Current numerical prototype

The repository now includes two small, auditable FEM verification paths.

### 1D Poisson

[
-\frac{d}{dx}\left(k\frac{du}{dx}\right)=f(x), \qquad x\in(0,L)
]

- Uniform 1D mesh.
- Linear P1 elements.
- Two-point Gauss quadrature.
- Strong endpoint Dirichlet conditions.
- CSR assembly and conjugate-gradient solution.

### 2D scalar diffusion

[
-\nabla\cdot(k\nabla u)=f \quad \text{in }\Omega
]

- Triangular P1 mesh and structured rectangular mesh generator.
- Positive six-point degree-four reference-triangle quadrature.
- Jacobian-based physical basis gradients.
- Positive scalar isotropic diffusivity, including spatially varying coefficients.
- Volumetric source integration.
- Strong nodal Dirichlet elimination.
- Linear harmonic patch benchmark.
- Manufactured quartic solution with a second-order nodal-error refinement test.

The 2D implementation is deliberately limited: it does not yet support Neumann/Robin conditions, pure Neumann systems, anisotropic diffusion, P2 elements, adaptive meshing, L2/H1 error norms, or a general DOF manager. These are planned engineering tasks, not existing features.

## Planned architecture

- **Application:** Qt desktop interface and project/model management.
- **Model layer:** geometry entities, materials, fields, parameters, units, boundary conditions, and studies.
- **Geometry and mesh:** CAD exchange, topology, mesh generation, validation, and refinement.
- **Numerical core:** function spaces, element definitions, quadrature, DOF mapping, assembly, sparse algebra, and solvers.
- **Physics:** thermal, electrostatic, structural, electromagnetic, fluid, acoustic, and semiconductor interfaces.
- **Multiphysics:** segregated and monolithic coupling.
- **Results:** datasets, derived quantities, visualization, and reporting.
- **Automation:** Python API and reproducible job execution.

## Repository layout

```text
NOVA-Multiphysics/
├── CMakeLists.txt
├── CMakePresets.json
├── README.md
├── .github/workflows/ci.yml
├── docs/
├── examples/
│   ├── poisson_1d/
│   └── diffusion_2d/
├── include/nova/
│   ├── linalg/
│   └── fem/
├── src/
│   ├── linalg/
│   └── fem/
└── tests/
```

## Build

Requirements:
- CMake 3.20 or newer
- A C++20 compiler (GCC 11+, Clang 14+, or MSVC 2022 recommended)

Configure and build:

```bash
cmake --preset default
cmake --build --preset default
ctest --preset default
```

If your generator does not support build presets, use:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Run the examples:

```bash
./build/nova_poisson_1d
./build/nova_diffusion_2d
```

On Windows, use the corresponding executable paths under the build directory (or the Release subdirectory for multi-configuration generators).

## Engineering principles

- Numerical correctness before feature breadth.
- Explicit units and documented sign conventions.
- Deterministic, reproducible tests.
- No silent fallback when a solver fails.
- Separate model definition from execution.
- Stable public interfaces and versioned project formats.
- Every physics feature must have analytical or independently validated benchmarks.
- Third-party dependencies must undergo license and security review before commercial distribution.

## Licensing

The repository's license is not yet finalized. Do not use this prototype as a basis for commercial engineering decisions until licensing, validation, and product-readiness reviews are complete.

## Disclaimer

NOVA is experimental research and software engineering work. It is not a certified engineering tool. Users must independently verify all results before relying on them.
