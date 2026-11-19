# Experimental C++ SDK install and consumer verification

NOVA's C++20 core can be installed as a relocatable CMake package. This is a **development SDK**, not a commercial, ABI-stable, or safety-certified release. The project's license is not finalized; do not infer redistribution rights from install support.

## Build and install

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNOVA_BUILD_TESTS=ON
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
cmake --install build --config Release --prefix ./nova-sdk
```

The install contains public `nova/*.hpp` headers, the `nova_core` compiled library and CMake package files. Examples and tests are not installed. On Windows multi-configuration generators use `--config Release` as shown.

## Use in a separate CMake project

```cmake
cmake_minimum_required(VERSION 3.20)
project(Consumer LANGUAGES CXX)
find_package(NOVA 0.8 CONFIG REQUIRED)
add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE NOVA::Core)
```

Configure with `-DCMAKE_PREFIX_PATH=/absolute/path/to/nova-sdk`, or set `NOVA_DIR` to its installed CMake package directory. Consumers inherit the C++20 requirement through the exported target.

For a reproducible compile/link/runtime test, configure and run `tests/package_consumer` in a *separate build directory* with that prefix. Hosted CI performs this check on Windows and Linux, not just a source-tree example build.

## Compatibility and release status

The package version is `0.8.0`, using same-minor package version compatibility while pre-1.0. There is no stable ABI guarantee. Follow the actual installed library's release criteria and license review before any commercial use. The main-branch SDK includes the verified steady 3D tetrahedral diffusion solver, Gmsh 2.2 importer, conduction flux postprocessing and VTK writer, together with the 1D and 2D APIs. Transient 3D heat studies are still developed separately and are not part of this package.
