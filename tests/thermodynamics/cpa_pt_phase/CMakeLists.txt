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

include("${CMAKE_CURRENT_LIST_DIR}/../../support/cmake/cpp_test_options.cmake")
foreach(target IN ITEMS mpmc_cpa_pt_phase_tests mpmc_cpa_pt_phase_headers)
    target_link_libraries(${target} PRIVATE mpmc::thermodynamics)
    mpmc_test_compile_options(${target} "${MPMC_CPA_PT_PHASE_ENABLE_SANITIZERS}"
        "CPA PT phase sanitizers require Linux GCC/Clang")
endforeach()

foreach(case IN ITEMS
        nonassociating_three_roots
        associating_helmholtz
        component_permutation
        selected_phase_fugacity
        selected_phase_derivatives
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
