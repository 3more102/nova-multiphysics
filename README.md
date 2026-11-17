# NOVA Multiphysics

**Engineering the Physics of Tomorrow.**

NOVA is an independent, C++20 multiphysics simulation platform under development. Its long-term objective is to provide an integrated environment for geometry, meshing, finite-element modeling, numerical solution, multiphysics coupling, and scientific visualization.

> **Development status:** Early engineering prototype. Not production-ready and not yet validated for engineering design decisions.

## Product direction

NOVA is being developed as a commercial, general-purpose computer-aided engineering (CAE) platform with an independently developed numerical core.

### Initial engineering priorities

1. Establish a reproducible C++20 build and test environment.
2. Implement numerical primitives and a small, auditable FEM foundation.
3. Verify each numerical method against analytical solutions and independent benchmarks.
4. Add physics interfaces through a stable module contract.
5. Build geometry, meshing, solver orchestration, visualization, and desktop workflows around the numerical core.

## Current prototype scope

The first milestone focuses on a one-dimensional Poisson boundary-value problem:

\[
-\frac{d}{dx}\left(k\frac{du}{dx}\right)=f(x), \qquad x\in(0,L)
\]

with prescribed Dirichlet boundary values. The benchmark uses constant \(k\), a uniform mesh, linear P1 elements, Gaussian quadrature, and a conjugate-gradient solver for the resulting symmetric positive-definite system.

This deliberately narrow benchmark is the first verification gate. It is not a claim of general-purpose FEM capability.

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
├── LICENSE
├── .gitignore
├── .github/workflows/ci.yml
├── docs/
│   ├── ARCHITECTURE.md
│   ├── NUMERICAL_CONVENTIONS.md
│   ├── ROADMAP.md
│   └── VERIFICATION_PLAN.md
├── examples/
│   └── poisson_1d/
├── include/nova/
│   ├── core/
│   ├── linalg/
│   ├── fem/
│   └── physics/
├── src/
│   ├── linalg/
│   ├── fem/
│   └── physics/
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
