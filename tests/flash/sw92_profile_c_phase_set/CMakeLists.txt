cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_profile_c_phase_set_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_profile_c_phase_set_tests
    sw92_profile_c_phase_set_test.cpp
    sw92_profile_c_phase_set_header.cpp)
target_link_libraries(mpmc_sw92_profile_c_phase_set_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_profile_c_phase_set_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_sw92_profile_c_phase_set_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_sw92_profile_c_phase_set_blocker_diagnostic
    blocker_diagnostic.cpp)
target_link_libraries(mpmc_sw92_profile_c_phase_set_blocker_diagnostic PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_profile_c_phase_set_blocker_diagnostic PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_sw92_profile_c_phase_set_blocker_diagnostic PROPERTIES
    CXX_EXTENSIONS OFF)
add_dependencies(mpmc_sw92_profile_c_phase_set_tests
    mpmc_sw92_profile_c_phase_set_blocker_diagnostic)

if(MSVC)
    target_compile_options(mpmc_sw92_profile_c_phase_set_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
    target_compile_options(mpmc_sw92_profile_c_phase_set_blocker_diagnostic PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_profile_c_phase_set_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    target_compile_options(mpmc_sw92_profile_c_phase_set_blocker_diagnostic PRIVATE
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
    target_compile_options(mpmc_sw92_profile_c_phase_set_blocker_diagnostic PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_profile_c_phase_set_blocker_diagnostic PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    publish_one_phase
    publish_two_phase
    publish_sample6_three_phase
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

add_test(NAME flash.sw92_profile_c_phase_set.blocker_diagnostic
         COMMAND mpmc_sw92_profile_c_phase_set_blocker_diagnostic)
set_tests_properties(flash.sw92_profile_c_phase_set.blocker_diagnostic PROPERTIES
    LABELS "flash;flash.sw92_profile_c_phase_set;diagnostic" TIMEOUT 300)
