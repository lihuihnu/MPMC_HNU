cmake_minimum_required(VERSION 3.21)
project(mpmc_flow_core_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/flow"
    flow-library)
add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad"
    ad-library)
enable_testing()

option(
    MPMC_FLOW_ENABLE_SANITIZERS
    "Enable ASan/UBSan for flow core tests"
    OFF)

add_executable(
    mpmc_flow_core_tests
    natural_variable_cell_state_test.cpp
    natural_variable_cell_state_header.cpp)
target_link_libraries(
    mpmc_flow_core_tests
    PRIVATE mpmc::flow)

add_executable(
    mpmc_flow_fugacity_tests
    fugacity_equilibrium_residual_test.cpp
    fugacity_equilibrium_residual_header.cpp)
target_link_libraries(
    mpmc_flow_fugacity_tests
    PRIVATE mpmc::flow mpmc::ad)

foreach(target IN ITEMS
        mpmc_flow_core_tests
        mpmc_flow_fugacity_tests)
    set_target_properties(
        ${target}
        PROPERTIES CXX_EXTENSIONS OFF)

    if(MSVC)
        target_compile_options(
            ${target}
            PRIVATE
                /W4
                /WX
                /permissive-
                /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(
            ${target}
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
                -Wconversion
                -Wshadow
                -Werror)
    endif()

    if(MPMC_FLOW_ENABLE_SANITIZERS)
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
            message(
                FATAL_ERROR
                "flow core sanitizers require Linux GCC/Clang")
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
endforeach()

set(core_cases
    layout
    valid_state
    invalid_state
    headers)
foreach(case IN LISTS core_cases)
    add_test(
        NAME flow.core.${case}
        COMMAND mpmc_flow_core_tests ${case})
    set_tests_properties(
        flow.core.${case}
        PROPERTIES
            LABELS "flow;natural_variables;contract"
            TIMEOUT 30)
endforeach()

set(fugacity_cases
    residual_values
    actual_phase_pressure
    differentiable_interface
    invalid_inputs
    headers)
foreach(case IN LISTS fugacity_cases)
    add_test(
        NAME flow.fugacity.${case}
        COMMAND mpmc_flow_fugacity_tests ${case})
    set_tests_properties(
        flow.fugacity.${case}
        PROPERTIES
            LABELS "flow;natural_variables;fugacity;thermodynamics_boundary"
            TIMEOUT 30)
endforeach()
