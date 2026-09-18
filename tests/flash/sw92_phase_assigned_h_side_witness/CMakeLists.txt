cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_phase_assigned_h_side_witness_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_phase_assigned_h_side_witness_tests
    sw92_phase_assigned_h_side_witness_test.cpp
    sw92_phase_assigned_h_side_witness_header.cpp)
target_link_libraries(mpmc_sw92_phase_assigned_h_side_witness_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_phase_assigned_h_side_witness_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/sw92")
set_target_properties(mpmc_sw92_phase_assigned_h_side_witness_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_phase_assigned_h_side_witness_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_phase_assigned_h_side_witness_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_PHASE_ASSIGNED_H_SIDE_WITNESS_ENABLE_SANITIZERS
       "Check SW92 Profile-C C2a1 with ASan/UBSan" OFF)
if(MPMC_SW92_PHASE_ASSIGNED_H_SIDE_WITNESS_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 C2a1 sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_phase_assigned_h_side_witness_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_phase_assigned_h_side_witness_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    co2_no_witness
    methane_brine_no_witness
    cross_family_separation
    synthetic_admissible_witness
    synthetic_role_guard
    resource_exhaustion
    root_failure
    source_guard
    permutation
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_phase_assigned_h_side_witness.${case}
             COMMAND mpmc_sw92_phase_assigned_h_side_witness_tests ${case})
    set_tests_properties(flash.sw92_phase_assigned_h_side_witness.${case} PROPERTIES
        LABELS "flash;flash.sw92_phase_assigned_h_side_witness;unit" TIMEOUT 120)
endforeach()
