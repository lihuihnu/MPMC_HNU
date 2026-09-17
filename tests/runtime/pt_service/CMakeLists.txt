cmake_minimum_required(VERSION 3.21)
project(mpmc_pt_service_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/runtime" runtime-library)
enable_testing()

option(MPMC_PT_SERVICE_ENABLE_SANITIZERS
       "Enable ASan/UBSan for PT service boundary tests" OFF)

add_executable(mpmc_pt_service_tests
    pt_service_test.cpp
    pt_service_contract_headers.cpp
    pt_service_headers.cpp)
target_link_libraries(mpmc_pt_service_tests PRIVATE mpmc::runtime)
set_target_properties(mpmc_pt_service_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pt_service_three_backend_tests
    three_backend_discovery_test.cpp)
target_link_libraries(mpmc_pt_service_three_backend_tests PRIVATE mpmc::runtime)
target_include_directories(mpmc_pt_service_three_backend_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_pt_service_three_backend_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pt_service_pr_sw_interchangeability_tests
    pr_sw_interchangeability_test.cpp)
target_link_libraries(mpmc_pt_service_pr_sw_interchangeability_tests PRIVATE mpmc::runtime)
target_include_directories(mpmc_pt_service_pr_sw_interchangeability_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_pt_service_pr_sw_interchangeability_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pt_service_sw92_co2_water_experimental_tests
    sw92_co2_water_experimental_regression_test.cpp)
target_link_libraries(mpmc_pt_service_sw92_co2_water_experimental_tests PRIVATE mpmc::runtime)
target_include_directories(mpmc_pt_service_sw92_co2_water_experimental_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_pt_service_sw92_co2_water_experimental_tests PROPERTIES CXX_EXTENSIONS OFF)

add_dependencies(mpmc_pt_service_tests
    mpmc_pt_service_three_backend_tests
    mpmc_pt_service_pr_sw_interchangeability_tests
    mpmc_pt_service_sw92_co2_water_experimental_tests)

set(all_targets
    mpmc_pt_service_tests
    mpmc_pt_service_three_backend_tests
    mpmc_pt_service_pr_sw_interchangeability_tests
    mpmc_pt_service_sw92_co2_water_experimental_tests)

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

if(MPMC_PT_SERVICE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PT service sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

set(cases
    discovery_inventory
    request_dispatch_provenance
    variable_phase_count
    scientific_outcomes
    request_errors
    backend_errors
    contract_guards
    configuration_guards
    headers)
foreach(case IN LISTS cases)
    add_test(NAME runtime.pt_service.${case}
             COMMAND mpmc_pt_service_tests ${case})
    set_tests_properties(runtime.pt_service.${case} PROPERTIES
        LABELS "runtime;service;pt_flash;contract" TIMEOUT 30)
endforeach()

add_test(NAME runtime.pt_service.three_backend_discovery
         COMMAND mpmc_pt_service_three_backend_tests)
set_tests_properties(runtime.pt_service.three_backend_discovery PROPERTIES
    LABELS "runtime;service;pt_flash;pr76;sw92;cpa;capability;conformance" TIMEOUT 30)

add_test(NAME runtime.pt_service.pr_sw_interchangeability
         COMMAND mpmc_pt_service_pr_sw_interchangeability_tests)
set_tests_properties(runtime.pt_service.pr_sw_interchangeability PROPERTIES
    LABELS "runtime;service;pt_flash;pr76;sw92;backend_switch;conformance" TIMEOUT 30)

add_test(NAME runtime.pt_service.sw92_co2_water_experimental
         COMMAND mpmc_pt_service_sw92_co2_water_experimental_tests)
set_tests_properties(runtime.pt_service.sw92_co2_water_experimental PROPERTIES
    LABELS "runtime;service;pt_flash;sw92;experimental;validation" TIMEOUT 30)
