cmake_minimum_required(VERSION 3.21)
project(mpmc_flow_discretization_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/flow/discretization"
    flow-discretization-library)
enable_testing()

option(
    MPMC_FLOW_DISCRETIZATION_ENABLE_SANITIZERS
    "Enable ASan/UBSan for flow-discretization tests"
    OFF)

set(
    MPMC_WELL_DISCRETIZATION_ENABLE_SANITIZERS
    "${MPMC_FLOW_DISCRETIZATION_ENABLE_SANITIZERS}")
add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../well/discretization"
    well-discretization-tests)

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

add_executable(
    mpmc_flow_discretization_normalized_spatial_tests
    normalized_component_face_contribution_test.cpp
    normalized_component_face_contribution_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_normalized_spatial_tests
    PRIVATE mpmc::flow_discretization)

add_executable(
    mpmc_flow_discretization_local_conservation_tests
    local_component_conservation_residual_test.cpp
    local_component_conservation_residual_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_local_conservation_tests
    PRIVATE mpmc::flow_discretization)

add_executable(
    mpmc_flow_discretization_multi_cell_conservation_tests
    owned_multi_cell_component_conservation_test.cpp
    owned_multi_cell_component_conservation_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_multi_cell_conservation_tests
    PRIVATE mpmc::flow_discretization)

add_executable(
    mpmc_flow_discretization_energy_tests
    energy_conservation_test.cpp
    energy_face_flux_header.cpp
    local_energy_conservation_residual_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_energy_tests
    PRIVATE mpmc::flow_discretization)

add_executable(
    mpmc_flow_discretization_cell_source_tests
    cell_source_test.cpp
    cell_source_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_cell_source_tests
    PRIVATE mpmc::flow_discretization)

foreach(target IN ITEMS
        mpmc_flow_discretization_tests
        mpmc_flow_discretization_component_flux_tests
        mpmc_flow_discretization_scatter_tests
        mpmc_flow_discretization_normalized_spatial_tests
        mpmc_flow_discretization_local_conservation_tests
        mpmc_flow_discretization_multi_cell_conservation_tests
        mpmc_flow_discretization_energy_tests
        mpmc_flow_discretization_cell_source_tests)
    set_target_properties(
        ${target}
        PROPERTIES CXX_EXTENSIONS OFF)
endforeach()

foreach(target IN ITEMS
        mpmc_flow_discretization_tests
        mpmc_flow_discretization_component_flux_tests
        mpmc_flow_discretization_scatter_tests
        mpmc_flow_discretization_normalized_spatial_tests
        mpmc_flow_discretization_local_conservation_tests
        mpmc_flow_discretization_multi_cell_conservation_tests
        mpmc_flow_discretization_energy_tests
        mpmc_flow_discretization_cell_source_tests)
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
            mpmc_flow_discretization_scatter_tests
            mpmc_flow_discretization_normalized_spatial_tests
            mpmc_flow_discretization_local_conservation_tests
            mpmc_flow_discretization_multi_cell_conservation_tests
            mpmc_flow_discretization_energy_tests
        mpmc_flow_discretization_cell_source_tests)
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

set(normalized_spatial_cases
    normalized_values_and_weighted_conservation
    normalized_jacobian
    jacobian_fresh_source_perturbation
    zero_primal_nonzero_jacobian
    invalid_inputs
    headers)

foreach(case IN LISTS normalized_spatial_cases)
    add_test(
        NAME flow_discretization.normalized_spatial.${case}
        COMMAND mpmc_flow_discretization_normalized_spatial_tests ${case})
    set_tests_properties(
        flow_discretization.normalized_spatial.${case}
        PROPERTIES
            LABELS "flow_discretization;component_flux;bulk_volume;normalized_spatial"
            TIMEOUT 30)
endforeach()

set(local_conservation_cases
    multi_face_residual_and_blocks
    jacobian_fresh_source_perturbation
    zero_primal_nonzero_jacobian
    invalid_inputs
    headers)

foreach(case IN LISTS local_conservation_cases)
    add_test(
        NAME flow_discretization.local_conservation.${case}
        COMMAND mpmc_flow_discretization_local_conservation_tests ${case})
    set_tests_properties(
        flow_discretization.local_conservation.${case}
        PROPERTIES
            LABELS "flow_discretization;component_conservation;local_residual;multi_face"
            TIMEOUT 30)
endforeach()

set(multi_cell_conservation_cases
    serial_patch_rows_and_conservation
    jacobian_mapping
    invalid_inputs
    headers)

foreach(case IN LISTS multi_cell_conservation_cases)
    add_test(
        NAME flow_discretization.multi_cell_conservation.${case}
        COMMAND mpmc_flow_discretization_multi_cell_conservation_tests ${case})
    set_tests_properties(
        flow_discretization.multi_cell_conservation.${case}
        PROPERTIES
            LABELS "flow_discretization;component_conservation;multi_cell;closed_patch"
            TIMEOUT 30)
endforeach()

set(energy_conservation_cases
    face_rate_and_fresh_perturbation
    normalized_local_and_closed_patch
    invalid_inputs
    headers)

foreach(case IN LISTS energy_conservation_cases)
    add_test(
        NAME flow_discretization.energy_conservation.${case}
        COMMAND mpmc_flow_discretization_energy_tests ${case})
    set_tests_properties(
        flow_discretization.energy_conservation.${case}
        PROPERTIES
            LABELS "flow_discretization;energy;enthalpy;conduction;conservation"
            TIMEOUT 30)
endforeach()

set(cell_source_cases
    normalized_sign_units_and_jacobian
    invalid_inputs
    headers)

foreach(case IN LISTS cell_source_cases)
    add_test(
        NAME flow_discretization.cell_source.${case}
        COMMAND mpmc_flow_discretization_cell_source_tests ${case})
    set_tests_properties(
        flow_discretization.cell_source.${case}
        PROPERTIES
            LABELS "flow_discretization;cell_source;component_conservation;energy"
            TIMEOUT 30)
endforeach()
