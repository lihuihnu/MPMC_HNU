cmake_minimum_required(VERSION 3.21)
project(mpmc_flow_discretization_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/flow_discretization"
    flow-discretization-library)
enable_testing()

option(
    MPMC_FLOW_DISCRETIZATION_ENABLE_SANITIZERS
    "Enable ASan/UBSan for flow-discretization tests"
    OFF)

add_executable(
    mpmc_flow_discretization_tests
    tpfa_phase_darcy_flux_test.cpp
    tpfa_phase_darcy_flux_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_tests
    PRIVATE mpmc::flow_discretization)
set_target_properties(
    mpmc_flow_discretization_tests
    PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(
        mpmc_flow_discretization_tests
        PRIVATE
            /W4
            /WX
            /permissive-
            /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(
        mpmc_flow_discretization_tests
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
            -Werror)
endif()

if(MPMC_FLOW_DISCRETIZATION_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(
            FATAL_ERROR
            "flow-discretization sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(
        mpmc_flow_discretization_tests
        PRIVATE
            -fsanitize=address,undefined
            -fno-sanitize-recover=all
            -fno-omit-frame-pointer)
    target_link_options(
        mpmc_flow_discretization_tests
        PRIVATE
            -fsanitize=address,undefined
            -fno-sanitize-recover=all)
endif()

set(cases
    materialized_phase_flux
    flux_jacobian_fresh_perturbation
    exact_zero_tie_branch
    gate_rejection
    headers)

foreach(case IN LISTS cases)
    add_test(
        NAME flow_discretization.tpfa_phase_flux.${case}
        COMMAND mpmc_flow_discretization_tests ${case})
    set_tests_properties(
        flow_discretization.tpfa_phase_flux.${case}
        PROPERTIES
            LABELS "flow_discretization;tpfa;phase_flux;darcy"
            TIMEOUT 30)
endforeach()
