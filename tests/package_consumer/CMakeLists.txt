cmake_minimum_required(VERSION 3.20)
project(NOVA_Package_Consumer LANGUAGES CXX)

find_package(NOVA 0.8 CONFIG REQUIRED)

add_executable(nova_package_consumer main.cpp)
target_link_libraries(nova_package_consumer PRIVATE NOVA::Core)

add_executable(nova_adaptive_refinement_consumer adaptive_refinement_consumer.cpp)
target_link_libraries(nova_adaptive_refinement_consumer PRIVATE NOVA::Core)

add_executable(nova_elasticity_consumer elasticity_consumer.cpp)
target_link_libraries(nova_elasticity_consumer PRIVATE NOVA::Core)

include(CTest)
if(BUILD_TESTING)
    add_test(NAME nova_installed_thermal_project
        COMMAND NOVA::Runner solve "${NOVA_PROJECT_EXAMPLE_DIR}/convecting_wall.nova")
    set_tests_properties(nova_installed_thermal_project PROPERTIES PASS_REGULAR_EXPRESSION "NOVA_RESULT 1")
    add_test(NAME nova_installed_package_consumer COMMAND nova_package_consumer)
    add_test(NAME nova_installed_elasticity_consumer COMMAND nova_elasticity_consumer)
    add_test(NAME nova_installed_adaptive_refinement_consumer COMMAND nova_adaptive_refinement_consumer)
endif()
