cmake_minimum_required(VERSION 3.21)
project(mpmc_flow_core_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/flow"
    flow-library)
if(NOT TARGET mpmc_ad)
    add_subdirectory(
        "${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad"
        ad-library)
endif()
if(NOT TARGET mpmc_flash)
    add_subdirectory(
        "${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash"
        flash-library)
endif()
enable_testing()

option(
    MPMC_FLOW_ENABLE_SANITIZERS
    "Enable ASan/UBSan for flow core tests"
    OFF)

add_executable(
    mpmc_flow_core_tests
    natural_variable_cell_state_test.cpp
    natural_variable_cell_state_header.cpp
    single_phase_natural_variable_test.cpp
    single_phase_natural_variable_header.cpp
    two_phase_natural_variable_test.cpp
    two_phase_natural_variable_header.cpp)
target_link_libraries(
    mpmc_flow_core_tests
    PRIVATE mpmc::flow)

add_executable(
    mpmc_flow_phase_transition_tests
    phase_set_transition_test.cpp
    phase_set_transition_header.cpp)
target_link_libraries(
    mpmc_flow_phase_transition_tests
    PRIVATE
        mpmc::flow
        mpmc::flash)

add_executable(
    mpmc_flow_fugacity_tests
    fugacity_equilibrium_residual_test.cpp
    fugacity_equilibrium_residual_header.cpp
    thermodynamics_fugacity_adapters_header.cpp)
target_link_libraries(
    mpmc_flow_fugacity_tests
    PRIVATE
        mpmc::flow_thermodynamics
        mpmc::ad)

add_executable(
    mpmc_flow_fugacity_linearization_tests
    fugacity_equilibrium_linearization_test.cpp
    fugacity_equilibrium_linearization_header.cpp)
target_link_libraries(
    mpmc_flow_fugacity_linearization_tests
    PRIVATE mpmc::flow)

add_executable(
    mpmc_flow_accumulation_tests
    component_accumulation_test.cpp
    component_accumulation_header.cpp)
target_link_libraries(
    mpmc_flow_accumulation_tests
    PRIVATE mpmc::flow)

add_executable(
    mpmc_flow_accumulation_time_tests
    component_accumulation_time_test.cpp
    component_accumulation_time_header.cpp)
target_link_libraries(
    mpmc_flow_accumulation_time_tests
    PRIVATE mpmc::flow)

add_executable(
    mpmc_flow_energy_accumulation_tests
    energy_accumulation_test.cpp
    energy_accumulation_header.cpp)
target_link_libraries(
    mpmc_flow_energy_accumulation_tests
    PRIVATE mpmc::flow)

add_executable(
    mpmc_flow_saturation_constitutive_tests
    saturation_constitutive_test.cpp
    saturation_constitutive_header.cpp)
target_link_libraries(
    mpmc_flow_saturation_constitutive_tests
    PRIVATE
        mpmc::flow
        mpmc::ad)

add_executable(
    mpmc_flow_phase_transport_tests
    phase_transport_test.cpp
    phase_transport_header.cpp)
target_link_libraries(
    mpmc_flow_phase_transport_tests
    PRIVATE
        mpmc::flow
        mpmc::ad)

add_executable(
    mpmc_flow_phase_potential_tests
    phase_potential_upwind_test.cpp
    phase_potential_upwind_header.cpp)
target_link_libraries(
    mpmc_flow_phase_potential_tests
    PRIVATE mpmc::flow)

foreach(target IN ITEMS
        mpmc_flow_core_tests
        mpmc_flow_phase_transition_tests
        mpmc_flow_fugacity_tests
        mpmc_flow_fugacity_linearization_tests
        mpmc_flow_accumulation_tests
        mpmc_flow_accumulation_time_tests
        mpmc_flow_energy_accumulation_tests
        mpmc_flow_saturation_constitutive_tests
        mpmc_flow_phase_transport_tests
        mpmc_flow_phase_potential_tests)
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
    composition_pivot
    valid_state
    invalid_state
    single_phase_reduction
    two_phase_reduction
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

add_test(
    NAME flow.phase_transition.contract
    COMMAND mpmc_flow_phase_transition_tests)
set_tests_properties(
    flow.phase_transition.contract
    PROPERTIES
        LABELS "flow;natural_variables;phase_transition;stability;conservation"
        TIMEOUT 30)

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

set(fugacity_linearization_cases
    carrier
    invalid_inputs
    headers)
foreach(case IN LISTS fugacity_linearization_cases)
    add_test(
        NAME flow.fugacity_linearization.${case}
        COMMAND mpmc_flow_fugacity_linearization_tests ${case})
    set_tests_properties(
        flow.fugacity_linearization.${case}
        PROPERTIES
            LABELS "flow;natural_variables;fugacity;linearization"
            TIMEOUT 30)
endforeach()

set(accumulation_cases
    snapshot_pair
    jacobian
    invalid_inputs
    headers)
foreach(case IN LISTS accumulation_cases)
    add_test(
        NAME flow.accumulation.${case}
        COMMAND mpmc_flow_accumulation_tests ${case})
    set_tests_properties(
        flow.accumulation.${case}
        PROPERTIES
            LABELS "flow;natural_variables;accumulation;conservation"
            TIMEOUT 30)
endforeach()

set(accumulation_time_cases
    backward_euler
    previous_is_frozen
    near_steady_roundoff
    invalid_inputs
    headers)
foreach(case IN LISTS accumulation_time_cases)
    add_test(
        NAME flow.accumulation_time.${case}
        COMMAND mpmc_flow_accumulation_time_tests ${case})
    set_tests_properties(
        flow.accumulation_time.${case}
        PROPERTIES
            LABELS "flow;natural_variables;accumulation;time_discretization;conservation"
            TIMEOUT 30)
endforeach()

set(energy_accumulation_cases
    accumulation_and_backward_euler
    jacobian_fresh_perturbation
    zero_residual_nonzero_jacobian
    invalid_inputs
    headers)
foreach(case IN LISTS energy_accumulation_cases)
    add_test(
        NAME flow.energy_accumulation.${case}
        COMMAND mpmc_flow_energy_accumulation_tests ${case})
    set_tests_properties(
        flow.energy_accumulation.${case}
        PROPERTIES
            LABELS "flow;natural_variables;energy;accumulation;conservation"
            TIMEOUT 30)
endforeach()

set(saturation_constitutive_cases
    primal
    differentiable_chart
    no_capillary
    invalid_inputs
    headers)
foreach(case IN LISTS saturation_constitutive_cases)
    add_test(
        NAME flow.constitutive.${case}
        COMMAND mpmc_flow_saturation_constitutive_tests ${case})
    set_tests_properties(
        flow.constitutive.${case}
        PROPERTIES
            LABELS "flow;natural_variables;relative_permeability;capillary_pressure;constitutive"
            TIMEOUT 30)
endforeach()

set(phase_transport_cases
    local_mobility
    no_capillary_pressure_gradient
    invalid_inputs
    headers)
foreach(case IN LISTS phase_transport_cases)
    add_test(
        NAME flow.transport.${case}
        COMMAND mpmc_flow_phase_transport_tests ${case})
    set_tests_properties(
        flow.transport.${case}
        PROPERTIES
            LABELS "flow;natural_variables;transport;mass_density;viscosity;mobility"
            TIMEOUT 30)
endforeach()

set(phase_potential_cases
    potential_and_upwind
    jacobian_fresh_perturbation
    orientation_reversal
    zero_gravity
    invalid_inputs
    headers)
foreach(case IN LISTS phase_potential_cases)
    add_test(
        NAME flow.potential.${case}
        COMMAND mpmc_flow_phase_potential_tests ${case})
    set_tests_properties(
        flow.potential.${case}
        PROPERTIES
            LABELS "flow;phase_potential;gravity;face_density;upwind"
            TIMEOUT 30)
endforeach()
