cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_asymmetric_orchestration_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_asymmetric_orchestration_tests
    sw92_asymmetric_orchestration_test.cpp
    sw92_asymmetric_orchestration_header.cpp)
target_link_libraries(mpmc_sw92_asymmetric_orchestration_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_asymmetric_orchestration_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_sw92_asymmetric_orchestration_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_asymmetric_orchestration_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_asymmetric_orchestration_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_ASYMMETRIC_ORCHESTRATION_ENABLE_SANITIZERS
       "Check SW92 Gate 3B.2 orchestration with ASan/UBSan" OFF)
if(MPMC_SW92_ASYMMETRIC_ORCHESTRATION_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 Gate 3B.2 sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_asymmetric_orchestration_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_asymmetric_orchestration_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    plan_and_same_family_alternatives
    full_candidate_selection
    family_neutral_attempt_budget
    stable_single_phase_candidate
    initial_family_tie_is_indeterminate
    slot_swap_deduplication
    distinct_candidate_gibbs_tie
    permutation
    resource_and_contract_failures
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_asymmetric_orchestration.${case}
             COMMAND mpmc_sw92_asymmetric_orchestration_tests ${case})
    set_tests_properties(flash.sw92_asymmetric_orchestration.${case} PROPERTIES
        LABELS "flash;flash.sw92_asymmetric_orchestration;unit" TIMEOUT 180)
endforeach()
