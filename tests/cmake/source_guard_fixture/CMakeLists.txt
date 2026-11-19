cmake_minimum_required(VERSION 3.20)
project(NOVASourceGuardFixture LANGUAGES CXX)

include("${CMAKE_CURRENT_LIST_DIR}/../../../cmake/NovaSourceGuard.cmake")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/fixture.cpp" "int fixture() { return 0; }\n")
add_library(guard_fixture STATIC "${CMAKE_CURRENT_BINARY_DIR}/fixture.cpp")
if(GUARD_DUPLICATE)
    target_sources(guard_fixture PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/fixture.cpp")
endif()
nova_assert_unique_target_sources(guard_fixture)
