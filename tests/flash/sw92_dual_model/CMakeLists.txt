cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_dual_model_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_dual_model_tests
    sw92_dual_model_test.cpp
    sw92_dual_model_header.cpp)
target_link_libraries(mpmc_sw92_dual_model_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_sw92_dual_model_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_dual_model_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_dual_model_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_DUAL_MODEL_ENABLE_SANITIZERS
       "Check SW92 dual-model orchestration with ASan/UBSan" OFF)
if(MPMC_SW92_DUAL_MODEL_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 dual-model sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_dual_model_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_dual_model_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    co2_dual_model
    methane_dual_model
    permutation
    family_failure_isolated
    multicomponent_label_guard
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_dual_model.${case}
             COMMAND mpmc_sw92_dual_model_tests ${case})
    set_tests_properties(flash.sw92_dual_model.${case} PROPERTIES
        LABELS "flash;flash.sw92_dual_model;unit" TIMEOUT 120)
endforeach()
