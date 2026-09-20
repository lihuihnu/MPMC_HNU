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

add_executable(
    mpmc_flow_discretization_component_flux_tests
    tpfa_component_molar_flux_test.cpp
    tpfa_component_molar_flux_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_component_flux_tests
    PRIVATE mpmc::flow_discretization)

add_executable(
    mpmc_flow_discretization_scatter_tests
    conservative_component_face_rate_scatter_test.cpp
    conservative_component_face_rate_scatter_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_scatter_tests
    PRIVATE mpmc::flow_discretization)
foreach(target IN ITEMS
        mpmc_flow_discretization_tests
        mpmc_flow_discretization_component_flux_tests
        mpmc_flow_discretization_scatter_tests)
    set_target_properties(
        ${target}
        PROPERTIES CXX_EXTENSIONS OFF)
endforeach()

foreach(target IN ITEMS
        mpmc_flow_discretization_tests
        mpmc_flow_discretization_component_flux_tests
        mpmc_flow_discretization_scatter_tests)
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
endforeach()

if(MPMC_FLOW_DISCRETIZATION_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(
            FATAL_ERROR
            "flow-discretization sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN ITEMS
            mpmc_flow_discretization_tests
            mpmc_flow_discretization_component_flux_tests
            mpmc_flow_discretization_scatter_tests)
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
    endforeach()
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

set(component_flux_cases
    upwind_phase_molar_content
    component_flux_and_closure
    jacobian_fresh_perturbation
    exact_zero_tie_uses_owner_content
    invalid_inputs
    headers)

foreach(case IN LISTS component_flux_cases)
    add_test(
        NAME flow_discretization.component_flux.${case}
        COMMAND mpmc_flow_discretization_component_flux_tests ${case})
    set_tests_properties(
        flow_discretization.component_flux.${case}
        PROPERTIES
            LABELS "flow_discretization;component_flux;molar;upwind"
            TIMEOUT 30)
endforeach()

set(scatter_cases
    conservative_values
    strict_jacobian_antisymmetry
    jacobian_fresh_source_perturbation
    zero_primal_nonzero_jacobian
    invalid_source
    headers)

foreach(case IN LISTS scatter_cases)
    add_test(
        NAME flow_discretization.scatter.${case}
        COMMAND mpmc_flow_discretization_scatter_tests ${case})
    set_tests_properties(
        flow_discretization.scatter.${case}
        PROPERTIES
            LABELS "flow_discretization;component_flux;conservation;scatter"
            TIMEOUT 30)
endforeach()
