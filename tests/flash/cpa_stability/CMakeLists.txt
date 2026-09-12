cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_stability_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_CPA_STABILITY_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA stability tests" OFF)

add_executable(mpmc_cpa_stability_tests
    cpa_stability_test.cpp
    cpa_stability_headers.cpp)
target_link_libraries(mpmc_cpa_stability_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_stability_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_cpa_stability_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_cpa_stability_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

if(MPMC_CPA_STABILITY_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA stability sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_cpa_stability_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_cpa_stability_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
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
