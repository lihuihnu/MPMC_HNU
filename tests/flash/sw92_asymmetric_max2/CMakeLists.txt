cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_asymmetric_max2_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_asymmetric_max2_tests
    sw92_asymmetric_max2_test.cpp
    sw92_asymmetric_max2_header.cpp)
target_link_libraries(mpmc_sw92_asymmetric_max2_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_asymmetric_max2_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/sw92")
set_target_properties(mpmc_sw92_asymmetric_max2_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_asymmetric_max2_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_asymmetric_max2_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_ASYMMETRIC_MAX2_ENABLE_SANITIZERS
       "Check SW92 asymmetric max2 final acceptance with ASan/UBSan" OFF)
if(MPMC_SW92_ASYMMETRIC_MAX2_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 asymmetric max2 sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_asymmetric_max2_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_asymmetric_max2_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    accepted_two_phase_reference
    final_common_tangent_and_required_starts
    accepted_single_phase
    selection_indeterminate_not_published
    final_indeterminate_not_published
    final_negative_witness_not_published
    pair_gibbs_above_feed_not_published
    publication_structural_guard
    component_permutation
    final_start_quota_preflight
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_asymmetric_max2.${case}
             COMMAND mpmc_sw92_asymmetric_max2_tests ${case})
    set_tests_properties(flash.sw92_asymmetric_max2.${case} PROPERTIES
        LABELS "flash;flash.sw92_asymmetric_max2;unit" TIMEOUT 180)
endforeach()
