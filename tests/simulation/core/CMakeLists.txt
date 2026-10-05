cmake_minimum_required(VERSION 3.21)
project(mpmc_simulation_core_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/simulation"
    simulation-library)
enable_testing()

option(
    MPMC_SIMULATION_ENABLE_SANITIZERS
    "Use AddressSanitizer and UndefinedBehaviorSanitizer for simulation core tests"
    OFF)

function(configure_simulation_test_target target)
    target_link_libraries(${target} PRIVATE mpmc::simulation)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(
            ${target}
            PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endif()
    if(MPMC_SIMULATION_ENABLE_SANITIZERS)
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
            message(FATAL_ERROR
                "Simulation core sanitizers require Linux GCC/Clang")
        endif()
        target_compile_options(
            ${target}
            PRIVATE -fsanitize=address,undefined -fno-sanitize-recover=all
                    -fno-omit-frame-pointer)
        target_link_options(
            ${target}
            PRIVATE -fsanitize=address,undefined -fno-sanitize-recover=all)
    endif()
endfunction()

add_executable(
    mpmc_simulation_timeline_tests
    simulation_timeline_test.cpp)
configure_simulation_test_target(mpmc_simulation_timeline_tests)

add_executable(
    mpmc_simulation_timeline_header_self_contained
    simulation_timeline_header_self_contained.cpp)
configure_simulation_test_target(mpmc_simulation_timeline_header_self_contained)

add_executable(
    mpmc_simulation_cursor_header_self_contained
    simulation_cursor_header_self_contained.cpp)
configure_simulation_test_target(mpmc_simulation_cursor_header_self_contained)

add_custom_target(
    mpmc_simulation_core_all
    DEPENDS
        mpmc_simulation_timeline_tests
        mpmc_simulation_timeline_header_self_contained
        mpmc_simulation_cursor_header_self_contained)

foreach(case IN ITEMS
    timeline_valid
    timeline_invalid
    cursor_initialization
    cap_composition
    accepted_progression
    boundary_skip_rejection
    terminal_state)
    add_test(
        NAME simulation.core.timeline.${case}
        COMMAND mpmc_simulation_timeline_tests ${case})
    set_tests_properties(
        simulation.core.timeline.${case}
        PROPERTIES
            LABELS "simulation;simulation.core;timeline;unit"
            TIMEOUT 30)
endforeach()

add_test(
    NAME simulation.core.timeline.header_self_contained
    COMMAND mpmc_simulation_timeline_header_self_contained)
set_tests_properties(
    simulation.core.timeline.header_self_contained
    PROPERTIES LABELS "simulation;simulation.core;timeline;header" TIMEOUT 30)

add_test(
    NAME simulation.core.cursor.header_self_contained
    COMMAND mpmc_simulation_cursor_header_self_contained)
set_tests_properties(
    simulation.core.cursor.header_self_contained
    PROPERTIES LABELS "simulation;simulation.core;cursor;header" TIMEOUT 30)
