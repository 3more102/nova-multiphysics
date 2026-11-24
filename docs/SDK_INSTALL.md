# Experimental C++ SDK install and consumer verification

NOVA's C++20 core can be installed as a relocatable CMake package. This is a **development SDK**, not a commercial, ABI-stable, or safety-certified release. The project's license is not finalized; do not infer redistribution rights from install support.

## Build and install

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNOVA_BUILD_TESTS=ON
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
cmake --install build --config Release --prefix ./nova-sdk
```

The install contains the public headers, compiled core library, CMake package, `nova-run` process and the native wall project. Test binaries and compiled numerical examples remain in the build tree. On Windows multi-configuration generators use `--config Release` as shown.

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

## Relocating a shared SDK

Configure with `-DBUILD_SHARED_LIBS=ON` to build a shared core. On Linux, the installed runner resolves the core relative to its executable location using an ELF origin-relative search path. On Windows the runtime DLL is installed beside `nova-run.exe`. Moving the entire SDK prefix preserves their relationship, including configurations that use relative, nondefault `CMAKE_INSTALL_BINDIR` and `CMAKE_INSTALL_LIBDIR`. Downstream Windows executables installed elsewhere need the SDK runtime directory on their DLL search path.

~~~sh
cmake -S . -B build-shared -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DNOVA_BUILD_EXAMPLES=OFF -DNOVA_BUILD_TESTS=OFF
cmake --build build-shared --config Release --parallel
python verification/verify_relocated_thermal_runner.py build-shared
~~~

The check installs to a temporary prefix containing spaces, moves it, runs an independent physical wall comparison and builds and tests a separate SDK consumer. On Linux it also resolves the loaded library to the moved installation.

## Compatibility and release status

The package version is `0.8.0`, using same-minor package version compatibility while pre-1.0. There is no stable ABI guarantee. Follow the actual installed library's release criteria and license review before any commercial use. The installed public headers include steady and transient 3D heat, Gmsh 2.2/4.1 import, conduction flux postprocessing, the VTK writer, the native stationary thermal project interface and the existing 1D and 2D APIs. See the [thermal project scope](THERMAL_PROJECT.md).
