cmake_minimum_required(VERSION 3.21)
project(mpmc_pr76_three_phase_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_PR76_THREE_PHASE_ENABLE_SANITIZERS
       "Enable ASan/UBSan for PR76 three-phase tests" OFF)

add_executable(mpmc_pr76_three_phase_tests
    pr76_three_phase_test.cpp
    pr76_three_phase_headers.cpp)
target_link_libraries(mpmc_pr76_three_phase_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_three_phase_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_three_phase_seeded_tests
    pr76_seeded_max3_test.cpp)
target_link_libraries(mpmc_pr76_three_phase_seeded_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_three_phase_seeded_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_three_phase_robustness_tests
    pr76_max3_robustness_test.cpp)
target_link_libraries(mpmc_pr76_three_phase_robustness_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_three_phase_robustness_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_three_phase_reference_test
    pr76_reference_candidate_test.cpp)
target_link_libraries(mpmc_pr76_three_phase_reference_test PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_three_phase_reference_test PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_sour_gas_three_phase_tests
    pr76_sour_gas_three_phase_test.cpp)
target_link_libraries(mpmc_pr76_sour_gas_three_phase_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_sour_gas_three_phase_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_heringer_2026_fixture_test
    pr76_heringer_2026_fixture_test.cpp)
target_link_libraries(mpmc_pr76_heringer_2026_fixture_test PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_heringer_2026_fixture_test PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_heringer_2026_boundary_test
    pr76_heringer_2026_boundary_test.cpp)
target_link_libraries(mpmc_pr76_heringer_2026_boundary_test PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_heringer_2026_boundary_test PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_strict_binary_kij_calibration_test
    pr76_strict_binary_kij_calibration_test.cpp)
target_link_libraries(mpmc_pr76_strict_binary_kij_calibration_test PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_strict_binary_kij_calibration_test PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_c10_kappa_mismatch_audit_test
    pr76_c10_kappa_mismatch_audit_test.cpp)
target_link_libraries(mpmc_pr76_c10_kappa_mismatch_audit_test PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_c10_kappa_mismatch_audit_test PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_c10_kij_structure_audit_test
    pr76_c10_kij_structure_audit_test.cpp)
target_link_libraries(mpmc_pr76_c10_kij_structure_audit_test PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_c10_kij_structure_audit_test PROPERTIES CXX_EXTENSIONS OFF)

add_dependencies(mpmc_pr76_three_phase_tests
    mpmc_pr76_three_phase_seeded_tests
    mpmc_pr76_three_phase_robustness_tests
    mpmc_pr76_three_phase_reference_test
    mpmc_pr76_sour_gas_three_phase_tests
    mpmc_pr76_heringer_2026_fixture_test
    mpmc_pr76_heringer_2026_boundary_test
    mpmc_pr76_strict_binary_kij_calibration_test
    mpmc_pr76_c10_kappa_mismatch_audit_test
    mpmc_pr76_c10_kij_structure_audit_test)

if(MSVC)
    foreach(target IN ITEMS
            mpmc_pr76_three_phase_tests
            mpmc_pr76_three_phase_seeded_tests
            mpmc_pr76_three_phase_robustness_tests
            mpmc_pr76_three_phase_reference_test
            mpmc_pr76_sour_gas_three_phase_tests
            mpmc_pr76_heringer_2026_fixture_test
            mpmc_pr76_heringer_2026_boundary_test
            mpmc_pr76_strict_binary_kij_calibration_test
            mpmc_pr76_c10_kappa_mismatch_audit_test
            mpmc_pr76_c10_kij_structure_audit_test)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    endforeach()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(target IN ITEMS
            mpmc_pr76_three_phase_tests
            mpmc_pr76_three_phase_seeded_tests
            mpmc_pr76_three_phase_robustness_tests
            mpmc_pr76_three_phase_reference_test
            mpmc_pr76_sour_gas_three_phase_tests
            mpmc_pr76_heringer_2026_fixture_test
            mpmc_pr76_heringer_2026_boundary_test
            mpmc_pr76_strict_binary_kij_calibration_test
            mpmc_pr76_c10_kappa_mismatch_audit_test
            mpmc_pr76_c10_kij_structure_audit_test)
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endforeach()
endif()

if(MPMC_PR76_THREE_PHASE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PR76 three-phase sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN ITEMS
            mpmc_pr76_three_phase_tests
            mpmc_pr76_three_phase_seeded_tests
            mpmc_pr76_three_phase_robustness_tests
            mpmc_pr76_three_phase_reference_test
            mpmc_pr76_sour_gas_three_phase_tests
            mpmc_pr76_heringer_2026_fixture_test
            mpmc_pr76_heringer_2026_boundary_test
            mpmc_pr76_strict_binary_kij_calibration_test
            mpmc_pr76_c10_kappa_mismatch_audit_test
            mpmc_pr76_c10_kij_structure_audit_test)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

set(cases
    generic_exact
    generic_disappearance
    pr76_fresh_two_phase_neighbor
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.pr76_three_phase.${case}
             COMMAND mpmc_pr76_three_phase_tests ${case})
    set_tests_properties(flash.pr76_three_phase.${case} PROPERTIES
        LABELS "flash;pr76;three_phase;max3" TIMEOUT 300)
endforeach()

foreach(case IN ITEMS
        seeded_max3
        seeded_component_permutation
        malformed_start
        start_requires_instability)
    add_test(NAME flash.pr76_three_phase.${case}
             COMMAND mpmc_pr76_three_phase_seeded_tests ${case})
    set_tests_properties(flash.pr76_three_phase.${case} PROPERTIES
        LABELS "flash;pr76;three_phase;max3;seeded_multistart" TIMEOUT 300)
endforeach()

foreach(case IN ITEMS
        automatic_cold_start
        continuation_budget
        continuation_3_to_3
        continuation_3_to_2)
    add_test(NAME flash.pr76_three_phase.${case}
             COMMAND mpmc_pr76_three_phase_robustness_tests ${case})
    set_tests_properties(flash.pr76_three_phase.${case} PROPERTIES
        LABELS "flash;pr76;three_phase;max3;robustness;continuation" TIMEOUT 300)
endforeach()

add_test(NAME flash.pr76_three_phase.reference_candidate
         COMMAND mpmc_pr76_three_phase_reference_test)
set_tests_properties(flash.pr76_three_phase.reference_candidate PROPERTIES
    LABELS "flash;pr76;three_phase;max3;reference_candidate" TIMEOUT 300)

add_test(NAME flash.pr76_three_phase.heringer_2026_fixture
         COMMAND mpmc_pr76_heringer_2026_fixture_test)
set_tests_properties(flash.pr76_three_phase.heringer_2026_fixture PROPERTIES
    LABELS "flash;pr76;three_phase;physical_mixture;literature_fixture" TIMEOUT 300)

add_test(NAME flash.pr76_three_phase.heringer_2026_3_to_2
         COMMAND mpmc_pr76_heringer_2026_boundary_test)
set_tests_properties(flash.pr76_three_phase.heringer_2026_3_to_2 PROPERTIES
    LABELS "flash;pr76;three_phase;max3;physical_mixture;literature_topology" TIMEOUT 300)

add_test(NAME flash.pr76_three_phase.strict_binary_kij_calibration
         COMMAND mpmc_pr76_strict_binary_kij_calibration_test)
set_tests_properties(flash.pr76_three_phase.strict_binary_kij_calibration PROPERTIES
    LABELS "flash;pr76;physical_mixture;experimental_binary;parameter_calibration" TIMEOUT 300)

add_test(NAME flash.pr76_three_phase.c10_kappa_mismatch_audit
         COMMAND mpmc_pr76_c10_kappa_mismatch_audit_test)
set_tests_properties(flash.pr76_three_phase.c10_kappa_mismatch_audit PROPERTIES
    LABELS "flash;pr76;physical_mixture;experimental_binary;model_convention_audit" TIMEOUT 300)

add_test(NAME flash.pr76_three_phase.c10_kij_structure_audit
         COMMAND mpmc_pr76_c10_kij_structure_audit_test)
set_tests_properties(flash.pr76_three_phase.c10_kij_structure_audit PROPERTIES
    LABELS "flash;pr76;physical_mixture;experimental_binary;mixing_rule_audit" TIMEOUT 300)

foreach(case IN ITEMS
        reference_candidate
        seeded_max3
        interior_reweighted
        component_permutation
        backend_publication)
    add_test(NAME flash.pr76_three_phase.sour_gas_${case}
             COMMAND mpmc_pr76_sour_gas_three_phase_tests ${case})
    set_tests_properties(flash.pr76_three_phase.sour_gas_${case} PROPERTIES
        LABELS "flash;pr76;three_phase;max3;physical_mixture;literature_reference" TIMEOUT 300)
endforeach()
