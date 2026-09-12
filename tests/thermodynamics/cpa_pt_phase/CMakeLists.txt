cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_pt_phase_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermodynamics-library)
enable_testing()

option(MPMC_CPA_PT_PHASE_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA PT phase tests" OFF)

add_executable(mpmc_cpa_pt_phase_tests
    cpa_pt_phase_test.cpp
    cpa_pt_phase_headers.cpp)
target_link_libraries(mpmc_cpa_pt_phase_tests PRIVATE mpmc::thermodynamics)
set_target_properties(mpmc_cpa_pt_phase_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_cpa_pt_phase_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_cpa_pt_phase_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

if(MPMC_CPA_PT_PHASE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA PT phase sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_cpa_pt_phase_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_cpa_pt_phase_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

foreach(case IN ITEMS
        nonassociating_three_roots
        associating_helmholtz
        component_permutation
        evaluation_budget)
    add_test(NAME thermodynamics.cpa_pt_phase.${case}
             COMMAND mpmc_cpa_pt_phase_tests ${case})
    set_tests_properties(thermodynamics.cpa_pt_phase.${case} PROPERTIES
        LABELS "thermodynamics;cpa;pt_phase;roots;fugacity" TIMEOUT 180)
endforeach()
