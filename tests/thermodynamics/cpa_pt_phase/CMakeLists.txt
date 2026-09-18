cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_pt_phase_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermodynamics-library)
enable_testing()

option(MPMC_CPA_PT_PHASE_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA PT phase tests" OFF)

add_executable(mpmc_cpa_pt_phase_tests
    cpa_pt_phase_test.cpp)
add_executable(mpmc_cpa_pt_phase_headers
    cpa_pt_phase_headers.cpp)

set(all_targets
    mpmc_cpa_pt_phase_tests
    mpmc_cpa_pt_phase_headers)
foreach(target IN LISTS all_targets)
    target_link_libraries(${target} PRIVATE mpmc::thermodynamics)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
endforeach()

if(MSVC)
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    endforeach()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endforeach()
endif()

if(MPMC_CPA_PT_PHASE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA PT phase sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
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

add_test(NAME thermodynamics.cpa_pt_phase.headers
         COMMAND mpmc_cpa_pt_phase_headers)
set_tests_properties(thermodynamics.cpa_pt_phase.headers PROPERTIES
    LABELS "thermodynamics;cpa;pt_phase;headers" TIMEOUT 30)
