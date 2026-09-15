cmake_minimum_required(VERSION 3.21)
project(mpmc_pt_flash_backend_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_PT_FLASH_BACKEND_ENABLE_SANITIZERS
       "Enable ASan/UBSan for PT flash backend tests" OFF)

add_executable(mpmc_pt_flash_backend_tests
    pt_flash_backend_test.cpp
    pt_flash_backend_headers.cpp)
target_link_libraries(mpmc_pt_flash_backend_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_pt_flash_backend_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_pt_flash_backend_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pt_flash_backend_pr76_max3_tests
    pr76_max3_seeded_backend_test.cpp)
target_link_libraries(mpmc_pt_flash_backend_pr76_max3_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pt_flash_backend_pr76_max3_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pt_flash_backend_pr76_near_dew_tests
    pr76_near_dew_backend_test.cpp)
target_link_libraries(mpmc_pt_flash_backend_pr76_near_dew_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pt_flash_backend_pr76_near_dew_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pt_flash_backend_cpa_max3_tests
    cpa_max3_backend_test.cpp)
target_link_libraries(mpmc_pt_flash_backend_cpa_max3_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pt_flash_backend_cpa_max3_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pt_flash_backend_freeze_tests
    three_backend_freeze_test.cpp)
target_link_libraries(mpmc_pt_flash_backend_freeze_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_pt_flash_backend_freeze_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_pt_flash_backend_freeze_tests PROPERTIES CXX_EXTENSIONS OFF)

add_dependencies(mpmc_pt_flash_backend_tests
    mpmc_pt_flash_backend_pr76_max3_tests
    mpmc_pt_flash_backend_pr76_near_dew_tests
    mpmc_pt_flash_backend_cpa_max3_tests
    mpmc_pt_flash_backend_freeze_tests)

set(all_targets
    mpmc_pt_flash_backend_tests
    mpmc_pt_flash_backend_pr76_max3_tests
    mpmc_pt_flash_backend_pr76_near_dew_tests
    mpmc_pt_flash_backend_cpa_max3_tests
    mpmc_pt_flash_backend_freeze_tests)

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

if(MPMC_PT_FLASH_BACKEND_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PT flash backend sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

set(cases
    sw92_direct_equivalence
    sw92_three_phase_capability
    sw92_boundary_transition_projection
    sw92_boundary_failure_projection
    transition_contract_guards
    contract_guards
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.pt_flash_backend.${case}
             COMMAND mpmc_pt_flash_backend_tests ${case})
    set_tests_properties(flash.pt_flash_backend.${case} PROPERTIES
        LABELS "flash;pt_flash_backend;phase_transition;boundary;capability;conformance" TIMEOUT 300)
endforeach()

add_test(NAME flash.pt_flash_backend.pr76_max3_equivalence
         COMMAND mpmc_pt_flash_backend_pr76_max3_tests)
set_tests_properties(flash.pt_flash_backend.pr76_max3_equivalence PROPERTIES
    LABELS "flash;pt_flash_backend;pr76;three_phase;max3;seeded_multistart;conformance" TIMEOUT 300)

add_test(NAME flash.pt_flash_backend.pr76_near_dew_classification
         COMMAND mpmc_pt_flash_backend_pr76_near_dew_tests)
set_tests_properties(flash.pt_flash_backend.pr76_near_dew_classification PROPERTIES
    LABELS "flash;pt_flash_backend;pr76;boundary;phase_disappearance;stability;conformance" TIMEOUT 300)

add_test(NAME flash.pt_flash_backend.cpa_max3_equivalence
         COMMAND mpmc_pt_flash_backend_cpa_max3_tests)
set_tests_properties(flash.pt_flash_backend.cpa_max3_equivalence PROPERTIES
    LABELS "flash;pt_flash_backend;cpa;three_phase;max3;seeded_multistart;conformance" TIMEOUT 300)

add_test(NAME flash.pt_flash_backend.three_backend_freeze
         COMMAND mpmc_pt_flash_backend_freeze_tests)
set_tests_properties(flash.pt_flash_backend.three_backend_freeze PROPERTIES
    LABELS "flash;pt_flash_backend;pr76;sw92;cpa;capability;contract_freeze;conformance" TIMEOUT 300)
