cmake_minimum_required(VERSION 3.20)

project(NOVA_Multiphysics
    VERSION 0.8.0
    DESCRIPTION "Independent multiphysics simulation platform"
    LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(cmake/NovaSourceGuard.cmake)

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
    src/linalg/incomplete_cholesky.cpp
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
    src/fem/transient_heat_3d.cpp
    src/fem/energy_balance_3d.cpp
    src/fem/conservative_face_heat_3d.cpp
    src/physics/electrostatics_3d.cpp
    src/io/gmsh22_tetrahedral.cpp
    src/io/gmsh_tetrahedral.cpp
    src/io/vtk_legacy_3d.cpp
    src/io/vtk_series.cpp
    src/io/thermal_checkpoint_3d.cpp
    src/fem/mesh_refinement_2d.cpp
    src/fem/mesh_refinement_3d.cpp
    src/fem/poisson_1d.cpp
    src/fem/diffusion_2d.cpp
    src/fem/diffusion_h1_2d.cpp
    src/fem/transient_heat_2d.cpp
    src/fem/adaptive_heat_2d.cpp
    src/fem/error_norms_2d.cpp
    src/fem/error_norms_3d.cpp
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

    add_executable(nova_transient_heat_3d examples/transient_heat_3d/main.cpp)
    target_link_libraries(nova_transient_heat_3d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_transient_heat_3d)

    add_executable(nova_gmsh_transient_heat_3d examples/gmsh_transient_heat_3d/main.cpp)
    target_link_libraries(nova_gmsh_transient_heat_3d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_gmsh_transient_heat_3d)

    add_executable(nova_electrostatics_3d examples/electrostatics_3d/main.cpp)
    target_link_libraries(nova_electrostatics_3d PRIVATE NOVA::Core)
    nova_enable_warnings(nova_electrostatics_3d)

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
        tests/test_incomplete_cholesky.cpp
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
        tests/test_conservative_face_heat_3d.cpp
        tests/test_electrostatics_3d.cpp
        tests/test_transient_heat_3d.cpp
        tests/test_transient_heat_3d_restart.cpp
        tests/test_gmsh22_tetrahedral.cpp
        tests/test_gmsh41_tetrahedral.cpp
        tests/test_vtk_legacy_3d.cpp
        tests/test_vtk_series.cpp
        tests/test_thermal_checkpoint_3d.cpp
        tests/test_transient_heat_3d_streaming.cpp
        tests/test_diffusion_h1_2d.cpp
        tests/test_transient_heat_2d.cpp
        tests/test_adaptive_heat_2d.cpp
        tests/test_error_norms_2d.cpp
        tests/test_error_norms_3d.cpp
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
    # Positive and negative configure tests verify the source-duplication guard.
    add_test(NAME nova_source_guard_accepts_unique COMMAND "${CMAKE_COMMAND}"
        -S "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/source_guard_fixture"
        -B "${CMAKE_CURRENT_BINARY_DIR}/source_guard_positive" -DGUARD_DUPLICATE=OFF)
    add_test(NAME nova_source_guard_rejects_duplicate COMMAND "${CMAKE_COMMAND}"
        -S "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/source_guard_fixture"
        -B "${CMAKE_CURRENT_BINARY_DIR}/source_guard_negative" -DGUARD_DUPLICATE=ON)
    set_tests_properties(nova_source_guard_rejects_duplicate PROPERTIES WILL_FAIL TRUE)
    if(NOVA_BUILD_EXAMPLES)
        add_test(NAME nova_diffusion_3d_smoke COMMAND nova_diffusion_3d)
        add_test(NAME nova_gmsh_diffusion_3d_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh" 7)
        add_test(NAME nova_gmsh_named_groups_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh"
            --material solid 1 --temperature fixed-temperature 0 --heat-source solid 20)
        add_test(NAME nova_gmsh41_named_groups_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron41.msh"
            --material solid 1 --temperature fixed-temperature 0)
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
        add_test(NAME nova_gmsh_mesh_quality_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh" --check-mesh)
        set_tests_properties(nova_gmsh_mesh_quality_smoke PROPERTIES
            PASS_REGULAR_EXPRESSION "Minimum mean ratio")
        add_test(NAME nova_gmsh_mesh_quality_gate_reject COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh"
            --check-mesh --min-mean-ratio 0.9)
        set_tests_properties(nova_gmsh_mesh_quality_gate_reject PROPERTIES WILL_FAIL TRUE)
        add_test(NAME nova_gmsh_mesh_quality_msh41_smoke COMMAND nova_gmsh_diffusion_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron41.msh" --check-mesh)
        set_tests_properties(nova_gmsh_mesh_quality_msh41_smoke PROPERTIES
            PASS_REGULAR_EXPRESSION "Minimum mean ratio")
        find_package(Python3 COMPONENTS Interpreter QUIET)
        if(Python3_Interpreter_FOUND)
            add_test(NAME nova_tagged_wall_verification COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/verify_tagged_conduction_3d.py"
                "$<TARGET_FILE:nova_gmsh_diffusion_3d>"
                "${CMAKE_CURRENT_SOURCE_DIR}/examples/gmsh_diffusion_3d/two_layer_wall.msh")
            add_test(NAME nova_gmsh_mesh_quality_analytical COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/verify_gmsh_quality_3d.py"
                "$<TARGET_FILE:nova_gmsh_diffusion_3d>"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh")
            add_test(NAME nova_gmsh_mesh_quality_analytical_msh41 COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/verify_gmsh_quality_3d.py"
                "$<TARGET_FILE:nova_gmsh_diffusion_3d>"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron41.msh")
        endif()
        add_test(NAME nova_transient_heat_3d_smoke COMMAND nova_transient_heat_3d)
        add_test(NAME nova_electrostatics_3d_smoke COMMAND nova_electrostatics_3d)
        set_tests_properties(nova_transient_heat_3d_smoke PROPERTIES
            PASS_REGULAR_EXPRESSION "NOVA 3D transient heat demonstration")
        if(Python3_Interpreter_FOUND)
            add_test(NAME nova_transient_heat_3d_pvd_verify COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/verify_transient_heat_3d_pvd.py"
                "$<TARGET_FILE:nova_transient_heat_3d>")
            add_test(NAME nova_gmsh_transient_heat_3d_verify
                COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/verify_gmsh_transient_heat_3d.py"
                "$<TARGET_FILE:nova_gmsh_transient_heat_3d>"
                "${CMAKE_CURRENT_SOURCE_DIR}/examples/gmsh_diffusion_3d/two_layer_wall.msh")
        endif()
        add_test(NAME nova_gmsh_transient_heat_3d_smoke
            COMMAND nova_gmsh_transient_heat_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/examples/gmsh_diffusion_3d/two_layer_wall.msh"
            --material 11 1 --material 22 3 --capacity 11 2 --capacity 22 2
            --initial 300 --source 4 --dt 0.1 --steps 5 --scheme cn)
        add_test(NAME nova_gmsh_transient_heat_3d_named_groups_smoke
            COMMAND nova_gmsh_transient_heat_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron.msh"
            --material solid 1 --capacity solid 1 --initial 300
            --temperature fixed-temperature 300 --dt 0.1 --steps 1)
        add_test(NAME nova_gmsh_transient_heat_3d_msh41_smoke
            COMMAND nova_gmsh_transient_heat_3d
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/unit_tetrahedron41.msh"
            --material solid 1 --capacity solid 1 --initial 320
            --temperature fixed-temperature 320 --dt 0.1 --steps 1)
        set_tests_properties(nova_gmsh_transient_heat_3d_msh41_smoke PROPERTIES
            PASS_REGULAR_EXPRESSION "min_temperature 320 max_temperature 320")
        add_test(NAME nova_transient_heat_smoke COMMAND nova_transient_heat_2d)
        add_test(NAME nova_transient_flux_smoke COMMAND nova_transient_flux_2d)
        add_test(NAME nova_adaptive_heat_smoke COMMAND nova_adaptive_heat_2d)
    endif()
endif()

# Fail configuration when a concurrent integration adds a source twice.
# A duplicate could otherwise hide changes in generated build files or CI.
nova_assert_unique_target_sources(nova_core)
if(NOVA_BUILD_TESTS)
    nova_assert_unique_target_sources(nova_tests)
endif()
if(NOVA_BUILD_EXAMPLES)
    foreach(example IN ITEMS
        nova_poisson_1d nova_diffusion_2d nova_diffusion_3d
        nova_gmsh_diffusion_3d nova_transient_heat_3d
        nova_gmsh_transient_heat_3d nova_electrostatics_3d nova_transient_heat_2d nova_transient_flux_2d
        nova_adaptive_heat_2d nova_nonlinear_diffusion_1d)
        nova_assert_unique_target_sources(${example})
    endforeach()
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
