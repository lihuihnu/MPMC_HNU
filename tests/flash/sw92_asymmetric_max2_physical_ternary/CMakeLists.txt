cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_asymmetric_max2_physical_ternary_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_asymmetric_max2_physical_ternary_tests
    sw92_asymmetric_max2_physical_ternary_test.cpp)
target_link_libraries(mpmc_sw92_asymmetric_max2_physical_ternary_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_asymmetric_max2_physical_ternary_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_sw92_asymmetric_max2_physical_ternary_tests PROPERTIES
    CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_asymmetric_max2_physical_ternary_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_asymmetric_max2_physical_ternary_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_MAX2_PHYSICAL_TERNARY_ENABLE_SANITIZERS
       "Check traceable SW92 max2 ternary regression with ASan/UBSan" OFF)
if(MPMC_SW92_MAX2_PHYSICAL_TERNARY_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 max2 ternary sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_asymmetric_max2_physical_ternary_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_asymmetric_max2_physical_ternary_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    traceable_parameter_snapshot
    ternary_max2_reference
    component_permutation
    missing_nonwater_pair_rejected)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_asymmetric_max2_physical_ternary.${case}
             COMMAND mpmc_sw92_asymmetric_max2_physical_ternary_tests ${case})
    set_tests_properties(flash.sw92_asymmetric_max2_physical_ternary.${case} PROPERTIES
        LABELS "flash;flash.sw92_asymmetric_max2_physical_ternary;validation" TIMEOUT 240)
endforeach()
