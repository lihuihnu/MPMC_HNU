cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_profile_c_sensitivity_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash/sensitivity" flash-sensitivity-library)
enable_testing()

option(MPMC_SW92_PROFILE_C_SENSITIVITY_ENABLE_SANITIZERS
       "Enable ASan/UBSan for SW92 Profile-C sensitivity tests" OFF)

add_executable(mpmc_sw92_profile_c_sensitivity_tests
    sw92_profile_c_sensitivity_test.cpp
    sw92_profile_c_sensitivity_header.cpp)
add_executable(mpmc_sw92_profile_c_sensitivity_diagnostic
    sw92_profile_c_sensitivity_diagnostic.cpp)
foreach(target IN ITEMS
        mpmc_sw92_profile_c_sensitivity_tests
        mpmc_sw92_profile_c_sensitivity_diagnostic)
    target_link_libraries(${target} PRIVATE mpmc::flash_sensitivity)
    target_include_directories(${target} PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endif()
endforeach()
add_dependencies(mpmc_sw92_profile_c_sensitivity_tests
    mpmc_sw92_profile_c_sensitivity_diagnostic)

if(MPMC_SW92_PROFILE_C_SENSITIVITY_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 Profile-C sensitivity sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN ITEMS
            mpmc_sw92_profile_c_sensitivity_tests
            mpmc_sw92_profile_c_sensitivity_diagnostic)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

set(cases
    phase_property_ad_crosscheck
    one_phase_fixed_set
    two_phase_fresh_resolve_crosscheck
    sample6_three_phase
    component_permutation_pt_columns
    h_slot_symmetry
    boundary_and_tamper_guards
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_profile_c_sensitivity.${case}
             COMMAND mpmc_sw92_profile_c_sensitivity_tests ${case})
    set_tests_properties(flash.sw92_profile_c_sensitivity.${case} PROPERTIES
        LABELS "flash;flash.sensitivity;sw92;ad" TIMEOUT 600)
endforeach()
add_test(NAME flash.sw92_profile_c_sensitivity.sample6_local_residual_diagnostic
         COMMAND mpmc_sw92_profile_c_sensitivity_diagnostic)
set_tests_properties(
    flash.sw92_profile_c_sensitivity.sample6_local_residual_diagnostic PROPERTIES
    LABELS "flash;flash.sensitivity;sw92;ad;diagnostic" TIMEOUT 600)
