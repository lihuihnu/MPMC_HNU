cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_asymmetric_max2_validation_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_asymmetric_max2_validation_tests
    sw92_asymmetric_max2_validation_test.cpp)
target_link_libraries(mpmc_sw92_asymmetric_max2_validation_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_asymmetric_max2_validation_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_sw92_asymmetric_max2_validation_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_asymmetric_max2_validation_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_asymmetric_max2_validation_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_MAX2_VALIDATION_ENABLE_SANITIZERS
       "Check SW92 max2 validation expansion with ASan/UBSan" OFF)
if(MPMC_SW92_MAX2_VALIDATION_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 max2 validation sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_asymmetric_max2_validation_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_asymmetric_max2_validation_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    ch4_fresh_max2
    ch4_brine_max2
    original_assignment_profile_separation
    runtime_ternary_structural
    phase_disappearance_boundary
    real_binary_family_envelope_audit
    synthetic_family_switch_guard
    root_envelope_nonsmooth)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_asymmetric_max2_validation.${case}
             COMMAND mpmc_sw92_asymmetric_max2_validation_tests ${case})
    set_tests_properties(flash.sw92_asymmetric_max2_validation.${case} PROPERTIES
        LABELS "flash;flash.sw92_asymmetric_max2_validation;validation" TIMEOUT 240)
endforeach()
