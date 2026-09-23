cmake_minimum_required(VERSION 3.21)
project(mpmc_well_discretization_tests LANGUAGES CXX)

if(NOT TARGET mpmc_well_discretization)
    add_subdirectory(
        "${CMAKE_CURRENT_LIST_DIR}/../../../modules/well/discretization"
        well-discretization-library)
endif()
enable_testing()

if(NOT DEFINED MPMC_WELL_DISCRETIZATION_ENABLE_SANITIZERS)
    option(
        MPMC_WELL_DISCRETIZATION_ENABLE_SANITIZERS
        "Enable ASan/UBSan for well-discretization tests"
        OFF)
endif()

function(configure_well_discretization_test_target target)
    target_link_libraries(
        ${target}
        PRIVATE mpmc::well_discretization)
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

    if(MPMC_WELL_DISCRETIZATION_ENABLE_SANITIZERS)
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
            message(
                FATAL_ERROR
                "Well-discretization sanitizers require Linux GCC/Clang")
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

add_executable(
    mpmc_well_discretization_tests
    hydraulic_conductance_test.cpp)
configure_well_discretization_test_target(
    mpmc_well_discretization_tests)

add_executable(
    mpmc_well_discretization_header_self_contained
    hydraulic_conductance_header_self_contained.cpp
    pressure_drawdown_rate_header_self_contained.cpp)
configure_well_discretization_test_target(
    mpmc_well_discretization_header_self_contained)

foreach(case IN ITEMS
    product_and_analytic_gradient
    phase_selection
    zero_mobility_nonzero_gradient
    pressure_drawdown_production_and_bhp_derivative
    pressure_drawdown_injection_sign_and_phase_selection
    pressure_drawdown_reservoir_product_rule
    pressure_drawdown_zero_drawdown
    pressure_drawdown_zero_mobility_nonzero_gradient
    pressure_drawdown_invalid_inputs
    invalid_inputs)
    add_test(
        NAME flow_discretization.well.hydraulic_conductance.${case}
        COMMAND mpmc_well_discretization_tests ${case})
    set_tests_properties(
        flow_discretization.well.hydraulic_conductance.${case}
        PROPERTIES
            LABELS "well;well_discretization;mobility;analytic_jacobian"
            TIMEOUT 30)
endforeach()

add_test(
    NAME flow_discretization.well.hydraulic_conductance.header_self_contained
    COMMAND mpmc_well_discretization_header_self_contained)
set_tests_properties(
    flow_discretization.well.hydraulic_conductance.header_self_contained
    PROPERTIES
        LABELS "well;well_discretization;header"
        TIMEOUT 30)
