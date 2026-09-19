cmake_minimum_required(VERSION 3.21)
project(mpmc_discretization_core_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/mesh"
    mesh-library)
add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/discretization"
    discretization-library)
enable_testing()

add_executable(mpmc_discretization_core_tests
    tpfa_test.cpp
    tpfa_half_connection_3d_header_self_contained.cpp
    tpfa_half_transmissibility_3d_header_self_contained.cpp
    tpfa_static_face_transmissibility_3d_header_self_contained.cpp
    tpfa_internal_face_transmissibility_snapshot_3d_header_self_contained.cpp)
target_link_libraries(
    mpmc_discretization_core_tests
    PRIVATE mpmc::discretization)
set_target_properties(
    mpmc_discretization_core_tests
    PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(
        mpmc_discretization_core_tests
        PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(
        mpmc_discretization_core_tests
        PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(
    MPMC_DISCRETIZATION_ENABLE_SANITIZERS
    "Use AddressSanitizer and UndefinedBehaviorSanitizer for discretization core tests"
    OFF)
if(MPMC_DISCRETIZATION_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR
            "Discretization core sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(
        mpmc_discretization_core_tests
        PRIVATE
            -fsanitize=address,undefined
            -fno-sanitize-recover=all
            -fno-omit-frame-pointer)
    target_link_options(
        mpmc_discretization_core_tests
        PRIVATE
            -fsanitize=address,undefined
            -fno-sanitize-recover=all)
endif()

foreach(case IN ITEMS
    tpfa_half_connection_3d
    tpfa_half_connection_3d_invalid
    tpfa_half_transmissibility_3d
    tpfa_half_transmissibility_3d_invalid
    tpfa_static_face_transmissibility_3d
    tpfa_static_face_transmissibility_3d_invalid
    tpfa_internal_face_transmissibility_snapshot_3d
    tpfa_internal_face_transmissibility_snapshot_3d_invalid)
    add_test(
        NAME discretization.core.${case}
        COMMAND mpmc_discretization_core_tests ${case})
    set_tests_properties(
        discretization.core.${case}
        PROPERTIES
            LABELS "discretization;discretization.core;tpfa;unit"
            TIMEOUT 30)
endforeach()
