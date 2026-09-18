cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_phase_assigned_boundary_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_phase_assigned_boundary_tests
    sw92_phase_assigned_boundary_test.cpp
    sw92_phase_assigned_boundary_header.cpp)
target_link_libraries(mpmc_sw92_phase_assigned_boundary_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_phase_assigned_boundary_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/sw92")
set_target_properties(mpmc_sw92_phase_assigned_boundary_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_phase_assigned_boundary_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_phase_assigned_boundary_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_PHASE_ASSIGNED_BOUNDARY_ENABLE_SANITIZERS
       "Check SW92 Profile-C boundary re-solve with ASan/UBSan" OFF)
if(MPMC_SW92_PHASE_ASSIGNED_BOUNDARY_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 boundary sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_phase_assigned_boundary_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_phase_assigned_boundary_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    physical_w_h_edge
    physical_no_w_edge
    synthetic_h_drop_guard
    synthetic_w_drop_guard
    source_chain_guard
    boundary_aware_interior
    boundary_aware_guard
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_phase_assigned_boundary.${case}
             COMMAND mpmc_sw92_phase_assigned_boundary_tests ${case})
    set_tests_properties(flash.sw92_phase_assigned_boundary.${case} PROPERTIES
        LABELS "flash;flash.sw92_phase_assigned_boundary;unit" TIMEOUT 240)
endforeach()
