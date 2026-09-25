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
    hydraulic_conductance_test.cpp
    component_molar_rate_test_support.cpp
    energy_rate_test_support.cpp
    cell_source_adapter_test_support.cpp)
configure_well_discretization_test_target(
    mpmc_well_discretization_tests)

add_executable(
    mpmc_well_discretization_header_self_contained
    hydraulic_conductance_header_self_contained.cpp
    pressure_drawdown_rate_header_self_contained.cpp)
configure_well_discretization_test_target(
    mpmc_well_discretization_header_self_contained)

add_executable(
    mpmc_well_component_molar_rate_header_self_contained
    component_molar_rate_header_self_contained.cpp)
configure_well_discretization_test_target(
    mpmc_well_component_molar_rate_header_self_contained)
add_dependencies(
    mpmc_well_discretization_header_self_contained
    mpmc_well_component_molar_rate_header_self_contained)

add_executable(
    mpmc_well_energy_rate_header_self_contained
    energy_rate_header_self_contained.cpp)
configure_well_discretization_test_target(
    mpmc_well_energy_rate_header_self_contained)
add_dependencies(
    mpmc_well_discretization_header_self_contained
    mpmc_well_energy_rate_header_self_contained)

add_executable(
    mpmc_well_cell_source_adapter_header_self_contained
    cell_source_adapter_header_self_contained.cpp)
configure_well_discretization_test_target(
    mpmc_well_cell_source_adapter_header_self_contained)
add_dependencies(
    mpmc_well_discretization_header_self_contained
    mpmc_well_cell_source_adapter_header_self_contained)

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
    component_molar_rate_mixed_phase_closure
    component_molar_rate_reservoir_product_rule_and_composition_pivot
    component_molar_rate_bhp_derivative_closure
    component_molar_rate_zero_phase_rate_preserves_derivative
    component_molar_rate_invalid_inputs
    energy_rate_production_and_injection_enthalpy_selection
    energy_rate_production_reservoir_enthalpy_product_rule
    energy_rate_injection_ignores_reservoir_enthalpy_derivative
    energy_rate_bhp_derivative_selected_enthalpy
    energy_rate_zero_tie_reservoir_enthalpy
    energy_rate_invalid_inputs
    cell_source_adapter_sign_and_reservoir_jacobian
    cell_source_adapter_bhp_sidecar
    cell_source_adapter_normalized_residual_sign
    cell_source_adapter_invalid_identity_and_inputs
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

add_test(
    NAME flow_discretization.well.component_molar_rate.header_self_contained
    COMMAND mpmc_well_component_molar_rate_header_self_contained)
set_tests_properties(
    flow_discretization.well.component_molar_rate.header_self_contained
    PROPERTIES
        LABELS "well;well_discretization;component_molar_rate;header"
        TIMEOUT 30)

add_test(
    NAME flow_discretization.well.energy_rate.header_self_contained
    COMMAND mpmc_well_energy_rate_header_self_contained)
set_tests_properties(
    flow_discretization.well.energy_rate.header_self_contained
    PROPERTIES
        LABELS "well;well_discretization;energy_rate;header"
        TIMEOUT 30)

add_test(
    NAME flow_discretization.well.cell_source_adapter.header_self_contained
    COMMAND mpmc_well_cell_source_adapter_header_self_contained)
set_tests_properties(
    flow_discretization.well.cell_source_adapter.header_self_contained
    PROPERTIES
        LABELS "well;well_discretization;cell_source_adapter;header"
        TIMEOUT 30)
