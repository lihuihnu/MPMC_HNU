cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_profile_c_phase_set_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
if(NOT TARGET mpmc::flow)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flow" flow-library)
endif()
enable_testing()

add_executable(mpmc_sw92_profile_c_phase_set_tests
    sw92_profile_c_phase_set_test.cpp
    sw92_profile_c_phase_set_header.cpp)
target_link_libraries(mpmc_sw92_profile_c_phase_set_tests PRIVATE mpmc::flash mpmc::flow_thermodynamics)
target_include_directories(mpmc_sw92_profile_c_phase_set_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/flow")
target_include_directories(mpmc_sw92_profile_c_phase_set_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/sw92")
set_target_properties(mpmc_sw92_profile_c_phase_set_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_sw92_profile_c_phase_set_provenance_tests
    sw92_profile_c_phase_set_provenance_test.cpp)
target_link_libraries(mpmc_sw92_profile_c_phase_set_provenance_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_profile_c_phase_set_provenance_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/sw92")
set_target_properties(mpmc_sw92_profile_c_phase_set_provenance_tests PROPERTIES
    CXX_EXTENSIONS OFF)
add_dependencies(mpmc_sw92_profile_c_phase_set_tests
    mpmc_sw92_profile_c_phase_set_provenance_tests)

if(MSVC)
    target_compile_options(mpmc_sw92_profile_c_phase_set_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
    target_compile_options(mpmc_sw92_profile_c_phase_set_provenance_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_profile_c_phase_set_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    target_compile_options(mpmc_sw92_profile_c_phase_set_provenance_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_PROFILE_C_PHASE_SET_ENABLE_SANITIZERS
       "Check SW92 Profile-C authoritative phase-set publication with ASan/UBSan" OFF)
if(MPMC_SW92_PROFILE_C_PHASE_SET_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 Profile-C phase-set sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_profile_c_phase_set_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_profile_c_phase_set_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
    target_compile_options(mpmc_sw92_profile_c_phase_set_provenance_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_profile_c_phase_set_provenance_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    publish_one_phase
    publish_two_phase
    publish_sample6_three_phase
    flow_fugacity_adapter
    fresh_disappearance_neighbors
    provenance_guard
    component_permutation
    h_slot_symmetry
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_profile_c_phase_set.${case}
             COMMAND mpmc_sw92_profile_c_phase_set_tests ${case})
    set_tests_properties(flash.sw92_profile_c_phase_set.${case} PROPERTIES
        LABELS "flash;flash.sw92_profile_c_phase_set;publication" TIMEOUT 300)
endforeach()

add_test(NAME flash.sw92_profile_c_phase_set.provenance_conventions
         COMMAND mpmc_sw92_profile_c_phase_set_provenance_tests)
set_tests_properties(flash.sw92_profile_c_phase_set.provenance_conventions PROPERTIES
    LABELS "flash;flash.sw92_profile_c_phase_set;publication;provenance" TIMEOUT 300)
