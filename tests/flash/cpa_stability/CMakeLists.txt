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

include("${CMAKE_CURRENT_LIST_DIR}/../../support/cmake/cpp_test_options.cmake")
foreach(target IN ITEMS mpmc_cpa_stability_tests mpmc_cpa_stability_headers)
    target_link_libraries(${target} PRIVATE mpmc::flash)
    mpmc_test_compile_options(${target} "${MPMC_CPA_STABILITY_ENABLE_SANITIZERS}"
        "CPA stability sanitizers require Linux GCC/Clang")
endforeach()

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
