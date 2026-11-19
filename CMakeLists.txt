cmake_minimum_required(VERSION 3.20)

project(NOVA_Multiphysics
    VERSION 0.8.0
    DESCRIPTION "Independent multiphysics simulation platform"
    LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

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
    src/fem/mesh3d.cpp
    src/fem/tetrahedron_h1.cpp
    src/fem/h1_space_3d.cpp
    src/fem/conductivity_3d.cpp
    src/fem/constrained_diffusion.cpp
    src/fem/diffusion_nullspace.cpp
    src/fem/diffusion_3d.cpp
    src/fem/p1_flux_3d.cpp
    src/fem/tagged_boundary_flux_3d.cpp
    src/fem/energy_balance_3d.cpp
    src/io/gmsh22_tetrahedral.cpp
    src/io/vtk_legacy_3d.cpp
    src/fem/mesh_refinement_2d.cpp
    src/fem/mesh_refinement_3d.cpp
    src/fem/poisson_1d.cpp
    src/fem/diffusion_2d.cpp
    src/fem/diffusion_h1_2d.cpp
    src/fem/transient_heat_2d.cpp
    src/fem/adaptive_heat_2d.cpp
    src/fem/error_norms_2d.cpp
    src/fem/nonlinear_diffusion_1d.cpp
    src/fem/h1_space_2d.cpp
)

add_library(NOVA::Core ALIAS nova_core)
target_compile_features(nova_core PUBLIC cxx_std_20)
target_include_directories(nova_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
set_target_properties(nova_core PROPERTIES
    CXX_EXTENSIONS OFF
    EXPORT_NAME Core
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

    add_executable(nova_diffusion_3d examples/diffusion_3d/main.cpp)
    target_link_libraries(nova_diffusion_3d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_diffusion_3d)

    add_executable(nova_gmsh_diffusion_3d examples/gmsh_diffusion_3d/main.cpp)
    target_link_libraries(nova_gmsh_diffusion_3d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_gmsh_diffusion_3d)

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
        tests/test_mesh_topology_3d.cpp
        tests/test_h1_space_3d.cpp
        tests/test_diffusion_3d.cpp
        tests/test_cell_conduction_3d.cpp
        tests/test_p1_flux_3d.cpp
        tests/test_tagged_boundary_flux_3d.cpp
        tests/test_energy_balance_3d.cpp
        tests/test_gmsh22_tetrahedral.cpp
        tests/test_vtk_legacy_3d.cpp
        tests/test_diffusion_h1_2d.cpp
        tests/test_transient_heat_2d.cpp
        tests/test_adaptive_heat_2d.cpp
        tests/test_error_norms_2d.cpp
        tests/test_mesh_topology_2d.cpp
        tests/test_mesh_refinement_2d.cpp
        tests/test_mesh_refinement_3d.cpp
        tests/test_nonlinear_diffusion_1d.cpp
        tests/test_h1_space_2d.cpp
    )
    target_include_directories(nova_tests PRIVATE tests)
    target_link_libraries(nova_tests PRIVATE NOVA::Core)
    nova_enable_warnings(nova_tests)
    add_test(NAME nova_unit_tests COMMAND nova_tests)
    if(NOVA_BUILD_EXAMPLES)
        add_test(NAME nova_diffusion_3d_smoke COMMAND nova_diffusion_3d)
        add_test(NAME nova_gmsh_diffusion_3d_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh" 7)
        add_test(NAME nova_gmsh_named_groups_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh"
            --material solid 1 --temperature fixed-temperature 0 --heat-source solid 20)
        add_test(NAME nova_gmsh_named_tensor_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh"
            --material-tensor solid 1 1 1 0 0 0 --temperature fixed-temperature 0)
        add_test(NAME nova_tagged_wall_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/gmsh_diffusion_3d/two_layer_wall.msh"
            --material 11 2 --material 22 6 --temperature 101 400
            --convection 102 12 300 --vtk "${CMAKE_CURRENT_BINARY_DIR}/two_layer_wall.vtk")
        add_test(NAME nova_tagged_heating_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/gmsh_diffusion_3d/two_layer_wall.msh"
            --material 11 2 --material 22 6 --temperature 101 400
            --convection 102 12 300 --heat-source 11 120 --heat-source 22 360)
        add_test(NAME nova_tagged_flux_report_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/gmsh_diffusion_3d/two_layer_wall.msh"
            --material 11 2 --material 22 6 --temperature 101 400
            --convection 102 12 300 --flux-report)
        set_tests_properties(nova_tagged_flux_report_smoke PROPERTIES
            PASS_REGULAR_EXPRESSION "Net outward conductive heat rate")
        add_test(NAME nova_energy_balance_3d_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/gmsh_diffusion_3d/two_layer_wall.msh"
            --material 11 2 --material 22 6 --temperature 101 400
            --convection 102 12 300 --balance-report)
        set_tests_properties(nova_energy_balance_3d_smoke PROPERTIES
            PASS_REGULAR_EXPRESSION "FEM signed energy balance defect")
        find_package(Python3 COMPONENTS Interpreter QUIET)
        if(Python3_Interpreter_FOUND)
            add_test(NAME nova_tagged_wall_verification COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/verify_tagged_conduction_3d.py"
                "$<TARGET_FILE:nova_gmsh_diffusion_3d>"
                "${CMAKE_CURRENT_SOURCE_DIR}/examples/gmsh_diffusion_3d/two_layer_wall.msh")
        endif()
        add_test(NAME nova_transient_heat_smoke COMMAND nova_transient_heat_2d)
        add_test(NAME nova_transient_flux_smoke COMMAND nova_transient_flux_2d)
        add_test(NAME nova_adaptive_heat_smoke COMMAND nova_adaptive_heat_2d)
    endif()
endif()

# Install a relocatable CMake package for downstream C++20 applications.
# This packages the experimental SDK, not a commercial release or license grant.
set(NOVA_INSTALL_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/NOVA")

install(TARGETS nova_core EXPORT NOVATargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/nova"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.hpp"
)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/NOVAConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/NOVAConfig.cmake"
    INSTALL_DESTINATION "${NOVA_INSTALL_CMAKE_DIR}"
)
# NOVA is pre-1.0: do not promise compatibility across minor versions.
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/NOVAConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMinorVersion
)
install(EXPORT NOVATargets
    FILE NOVATargets.cmake
    NAMESPACE NOVA::
    DESTINATION "${NOVA_INSTALL_CMAKE_DIR}"
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/NOVAConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/NOVAConfigVersion.cmake"
    DESTINATION "${NOVA_INSTALL_CMAKE_DIR}"
)
