cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_thermodynamic_closure_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/physics" physics-library)
enable_testing()

option(MPMC_SW92_THERMODYNAMIC_CLOSURE_ENABLE_SANITIZERS
       "Enable ASan/UBSan for SW92 thermodynamic-closure tests" OFF)

add_executable(mpmc_sw92_thermodynamic_closure_tests
    sw92_thermodynamic_closure_test.cpp
    sw92_thermodynamic_closure_header.cpp)
target_link_libraries(mpmc_sw92_thermodynamic_closure_tests PRIVATE mpmc::physics)
target_include_directories(mpmc_sw92_thermodynamic_closure_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/sw92")
set_target_properties(mpmc_sw92_thermodynamic_closure_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_thermodynamic_closure_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_thermodynamic_closure_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

if(MPMC_SW92_THERMODYNAMIC_CLOSURE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 thermodynamic-closure sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_thermodynamic_closure_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_thermodynamic_closure_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    one_phase_reproduces_missing_z
    two_phase_preserves_authoritative_properties
    sample6_three_phase_primal
    linearization_failure_policy
    rejection_and_model_guard
    component_permutation
    headers)
foreach(case IN LISTS cases)
    add_test(NAME physics.sw92_thermodynamic_closure.${case}
             COMMAND mpmc_sw92_thermodynamic_closure_tests ${case})
    set_tests_properties(physics.sw92_thermodynamic_closure.${case} PROPERTIES
        LABELS "physics;thermodynamic_closure;sw92;linearization" TIMEOUT 300)
endforeach()
