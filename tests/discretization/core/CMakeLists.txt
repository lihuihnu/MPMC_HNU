cmake_minimum_required(VERSION 3.21)
project(mpmc_discretization_core_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/mesh"
    mesh-library)
add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/discretization"
    discretization-library)
enable_testing()

option(
    MPMC_DISCRETIZATION_ENABLE_SANITIZERS
    "Use AddressSanitizer and UndefinedBehaviorSanitizer for discretization core tests"
    OFF)

set(
    MPMC_WELL_ENABLE_SANITIZERS
    "${MPMC_DISCRETIZATION_ENABLE_SANITIZERS}")
add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../well/core"
    well-core-tests)

function(configure_discretization_test_target target)
    target_link_libraries(
        ${target}
        PRIVATE mpmc::discretization)
    set_target_properties(
        ${target}
        PROPERTIES CXX_EXTENSIONS OFF)

    if(MSVC)
        target_compile_options(
            ${target}
            PRIVATE /W4 /WX /permissive- /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(
            ${target}
            PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endif()

    if(MPMC_DISCRETIZATION_ENABLE_SANITIZERS)
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
            message(FATAL_ERROR
                "Discretization core sanitizers require Linux GCC/Clang")
        endif()
        target_compile_options(
            ${target}
            PRIVATE
                -fsanitize=address,undefined
                -fno-sanitize-recover=all
                -fno-omit-frame-pointer)
        target_link_options(
            ${target}
            PRIVATE
                -fsanitize=address,undefined
                -fno-sanitize-recover=all)
    endif()
endfunction()

add_executable(mpmc_discretization_core_tests
    tpfa_test.cpp
    tpfa_half_connection_3d_header_self_contained.cpp
    tpfa_half_transmissibility_3d_header_self_contained.cpp
    tpfa_static_face_transmissibility_3d_header_self_contained.cpp
    tpfa_internal_face_transmissibility_snapshot_3d_header_self_contained.cpp)
configure_discretization_test_target(
    mpmc_discretization_core_tests)

add_executable(mpmc_discretization_admissibility_tests
    admissibility_test.cpp
    transmissibility_admissibility_3d_header_self_contained.cpp)
configure_discretization_test_target(
    mpmc_discretization_admissibility_tests)

add_custom_target(
    mpmc_discretization_core_all
    DEPENDS
        mpmc_discretization_core_tests
        mpmc_discretization_admissibility_tests
        mpmc_well_core_all)

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

foreach(case IN ITEMS
    geometry_admissibility_orthogonal
    geometry_admissibility_skewed
    k_orthogonality_axis_aligned
    k_orthogonality_skewed
    k_orthogonality_invalid)
    add_test(
        NAME discretization.core.${case}
        COMMAND mpmc_discretization_admissibility_tests ${case})
    set_tests_properties(
        discretization.core.${case}
        PROPERTIES
            LABELS "discretization;discretization.core;admissibility;unit"
            TIMEOUT 30)
endforeach()

