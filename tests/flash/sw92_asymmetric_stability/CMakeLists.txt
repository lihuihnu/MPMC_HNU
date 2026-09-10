cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_asymmetric_stability_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_asymmetric_stability_tests
    sw92_asymmetric_stability_test.cpp
    sw92_asymmetric_stability_header.cpp)
target_link_libraries(mpmc_sw92_asymmetric_stability_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_sw92_asymmetric_stability_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_asymmetric_stability_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_asymmetric_stability_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_ASYMMETRIC_STABILITY_ENABLE_SANITIZERS
       "Check SW92 asymmetric stability with ASan/UBSan" OFF)
if(MPMC_SW92_ASYMMETRIC_STABILITY_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 asymmetric stability sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_asymmetric_stability_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_asymmetric_stability_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    feed_reference_and_common_tangent
    negative_witness_and_counts
    pure_vertex_active_gauge
    family_tie_guard
    permutation
    family_failure_isolated
    same_family_nonsmooth_propagates
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_asymmetric_stability.${case}
             COMMAND mpmc_sw92_asymmetric_stability_tests ${case})
    set_tests_properties(flash.sw92_asymmetric_stability.${case} PROPERTIES
        LABELS "flash;flash.sw92_asymmetric_stability;unit" TIMEOUT 120)
endforeach()
