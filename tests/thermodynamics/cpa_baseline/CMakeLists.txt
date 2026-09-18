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

include("${CMAKE_CURRENT_LIST_DIR}/../../support/cmake/cpp_test_options.cmake")
foreach(target IN ITEMS mpmc_cpa_baseline_tests mpmc_cpa_baseline_headers)
    target_link_libraries(${target} PRIVATE mpmc::thermodynamics)
    mpmc_test_compile_options(${target} "${MPMC_CPA_BASELINE_ENABLE_SANITIZERS}"
        "CPA baseline sanitizers require Linux GCC/Clang")
endforeach()

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
