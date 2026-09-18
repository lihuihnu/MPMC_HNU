cmake_minimum_required(VERSION 3.21)
project(mpmc_component_inventory_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/physics" physics-library)
enable_testing()

option(MPMC_COMPONENT_INVENTORY_ENABLE_SANITIZERS
       "Enable ASan/UBSan for component-inventory tests" OFF)

add_executable(mpmc_component_inventory_tests
    component_inventory_test.cpp
    component_inventory_header.cpp)
target_link_libraries(mpmc_component_inventory_tests PRIVATE mpmc::physics)
target_include_directories(mpmc_component_inventory_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/sw92")
set_target_properties(mpmc_component_inventory_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_component_inventory_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_component_inventory_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

if(MPMC_COMPONENT_INVENTORY_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "component-inventory sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_component_inventory_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_component_inventory_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    synthetic_phase_permutation
    pr76_fixed_vle
    sw92_one_and_three_phase
    unavailable_policy
    zero_d_pr76
    zero_d_sw92
    headers)
foreach(case IN LISTS cases)
    add_test(NAME physics.component_inventory.${case}
             COMMAND mpmc_component_inventory_tests ${case})
    set_tests_properties(physics.component_inventory.${case} PROPERTIES
        LABELS "physics;component_inventory;linearization;zero_d" TIMEOUT 300)
endforeach()
