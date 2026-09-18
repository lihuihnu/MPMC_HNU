cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_baseline_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermodynamics-library)
enable_testing()

option(MPMC_CPA_BASELINE_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA baseline tests" OFF)

add_executable(mpmc_cpa_baseline_tests
    cpa_baseline_test.cpp)
add_executable(mpmc_cpa_baseline_headers
    cpa_headers.cpp)

set(all_targets
    mpmc_cpa_baseline_tests
    mpmc_cpa_baseline_headers)
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

if(MPMC_CPA_BASELINE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA baseline sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

foreach(case IN ITEMS
        parameter_reordering
        runtime_subset
        provenance_retention
        analytic_one_site
        association_pressure
        nonassociating_srk
        association_permutation
        missing_binary)
    add_test(NAME thermodynamics.cpa_baseline.${case}
             COMMAND mpmc_cpa_baseline_tests ${case})
    set_tests_properties(thermodynamics.cpa_baseline.${case} PROPERTIES
        LABELS "thermodynamics;cpa;association;baseline" TIMEOUT 120)
endforeach()

add_test(NAME thermodynamics.cpa_baseline.headers
         COMMAND mpmc_cpa_baseline_headers)
set_tests_properties(thermodynamics.cpa_baseline.headers PROPERTIES
    LABELS "thermodynamics;cpa;association;baseline;headers" TIMEOUT 30)
