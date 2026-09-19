cmake_minimum_required(VERSION 3.21)
project(mpmc_flow_core_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/flow"
    flow-library)
enable_testing()

option(
    MPMC_FLOW_ENABLE_SANITIZERS
    "Enable ASan/UBSan for flow core tests"
    OFF)

add_executable(
    mpmc_flow_core_tests
    natural_variable_cell_state_test.cpp
    natural_variable_cell_state_header.cpp)
target_link_libraries(
    mpmc_flow_core_tests
    PRIVATE mpmc::flow)
set_target_properties(
    mpmc_flow_core_tests
    PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(
        mpmc_flow_core_tests
        PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(
        mpmc_flow_core_tests
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
            -Werror)
endif()

if(MPMC_FLOW_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(
            FATAL_ERROR
            "flow core sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(
        mpmc_flow_core_tests
        PRIVATE
            -fsanitize=address,undefined
            -fno-sanitize-recover=all
            -fno-omit-frame-pointer)
    target_link_options(
        mpmc_flow_core_tests
        PRIVATE
            -fsanitize=address,undefined
            -fno-sanitize-recover=all)
endif()

set(cases
    layout
    valid_state
    invalid_state
    headers)
foreach(case IN LISTS cases)
    add_test(
        NAME flow.core.${case}
        COMMAND mpmc_flow_core_tests ${case})
    set_tests_properties(
        flow.core.${case}
        PROPERTIES
            LABELS "flow;natural_variables;contract"
            TIMEOUT 30)
endforeach()
