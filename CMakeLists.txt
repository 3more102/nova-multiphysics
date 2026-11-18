cmake_minimum_required(VERSION 3.20)

project(NOVA_Multiphysics
    VERSION 0.8.0
    DESCRIPTION "Independent multiphysics simulation platform"
    LANGUAGES CXX
)

option(NOVA_BUILD_TESTS "Build NOVA automated tests" ON)
option(NOVA_BUILD_EXAMPLES "Build NOVA examples" ON)

function(nova_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow
        )
    endif()
endfunction()

add_library(nova_core
    src/linalg/csr_matrix.cpp
    src/linalg/conjugate_gradient.cpp
    src/linalg/bicgstab.cpp
    src/linalg/gmres.cpp
    src/nonlinear/newton_solver.cpp
    src/fem/mesh2d.cpp
    src/fem/poisson_1d.cpp
    src/fem/diffusion_2d.cpp
    src/fem/diffusion_h1_2d.cpp
    src/fem/transient_heat_2d.cpp
    src/fem/adaptive_heat_2d.cpp
    src/fem/error_norms_2d.cpp
    src/fem/nonlinear_diffusion_1d.cpp
    src/fem/h1_space_2d.cpp
    src/io/vtk_legacy_2d.cpp
)

add_library(NOVA::Core ALIAS nova_core)
target_compile_features(nova_core PUBLIC cxx_std_20)
target_include_directories(nova_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
set_target_properties(nova_core PROPERTIES
    CXX_EXTENSIONS OFF
    POSITION_INDEPENDENT_CODE ON
)
nova_enable_warnings(nova_core)

if(NOVA_BUILD_EXAMPLES)
    add_executable(nova_poisson_1d examples/poisson_1d/main.cpp)
    target_link_libraries(nova_poisson_1d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_poisson_1d)

    add_executable(nova_diffusion_2d examples/diffusion_2d/main.cpp)
    target_link_libraries(nova_diffusion_2d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_diffusion_2d)

    add_executable(nova_transient_heat_2d examples/transient_heat_2d/main.cpp)
    target_link_libraries(nova_transient_heat_2d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_transient_heat_2d)

    add_executable(nova_transient_flux_2d examples/transient_flux_2d/main.cpp)
    target_link_libraries(nova_transient_flux_2d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_transient_flux_2d)

    add_executable(nova_adaptive_heat_2d examples/adaptive_heat_2d/main.cpp)
    target_link_libraries(nova_adaptive_heat_2d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_adaptive_heat_2d)

    add_executable(nova_nonlinear_diffusion_1d examples/nonlinear_diffusion_1d/main.cpp)
    target_link_libraries(nova_nonlinear_diffusion_1d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_nonlinear_diffusion_1d)
endif()

if(NOVA_BUILD_TESTS)
    enable_testing()
    add_executable(nova_tests
        tests/test_main.cpp
        tests/test_csr_matrix.cpp
        tests/test_conjugate_gradient.cpp
        tests/test_bicgstab.cpp
        tests/test_gmres.cpp
        tests/test_newton_solver.cpp
        tests/test_poisson_1d.cpp
        tests/test_diffusion_2d.cpp
        tests/test_diffusion_h1_2d.cpp
        tests/test_transient_heat_2d.cpp
        tests/test_adaptive_heat_2d.cpp
        tests/test_error_norms_2d.cpp
        tests/test_mesh_topology_2d.cpp
        tests/test_nonlinear_diffusion_1d.cpp
        tests/test_h1_space_2d.cpp
        tests/test_vtk_legacy_2d.cpp
    )
    target_include_directories(nova_tests PRIVATE tests)
    target_link_libraries(nova_tests PRIVATE NOVA::Core)
    nova_enable_warnings(nova_tests)
    add_test(NAME nova_unit_tests COMMAND nova_tests)
    if(NOVA_BUILD_EXAMPLES)
        add_test(NAME nova_transient_heat_smoke COMMAND nova_transient_heat_2d)
        add_test(NAME nova_transient_flux_smoke COMMAND nova_transient_flux_2d)
        add_test(NAME nova_adaptive_heat_smoke COMMAND nova_adaptive_heat_2d)
    endif()
endif()
