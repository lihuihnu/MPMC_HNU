cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_stability_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_CPA_STABILITY_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA stability tests" OFF)

add_executable(mpmc_cpa_stability_tests
    cpa_stability_test.cpp)
add_executable(mpmc_cpa_stability_headers
    cpa_stability_headers.cpp)

set(all_targets
    mpmc_cpa_stability_tests
    mpmc_cpa_stability_headers)
foreach(target IN LISTS all_targets)
    target_link_libraries(${target} PRIVATE mpmc::flash)
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

if(MPMC_CPA_STABILITY_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA stability sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

foreach(case IN ITEMS
        minimum_gibbs_root
        stable_reference
        negative_tpd
        component_permutation
        root_budget_failure)
    add_test(NAME flash.cpa_stability.${case}
             COMMAND mpmc_cpa_stability_tests ${case})
    set_tests_properties(flash.cpa_stability.${case} PROPERTIES
        LABELS "flash;cpa;stability;tpd" TIMEOUT 300)
endforeach()

add_test(NAME flash.cpa_stability.headers
         COMMAND mpmc_cpa_stability_headers)
set_tests_properties(flash.cpa_stability.headers PROPERTIES
    LABELS "flash;cpa;stability;headers" TIMEOUT 30)
