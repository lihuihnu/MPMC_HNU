cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_three_phase_experimental_line_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_three_phase_experimental_line_tests
    sw92_three_phase_experimental_line_test.cpp)
target_link_libraries(mpmc_sw92_three_phase_experimental_line_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_sw92_three_phase_experimental_line_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../support/sw92")
set_target_properties(mpmc_sw92_three_phase_experimental_line_tests PROPERTIES
    CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_three_phase_experimental_line_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_three_phase_experimental_line_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_THREE_PHASE_EXPERIMENTAL_LINE_ENABLE_SANITIZERS
       "Check SW92 experimental three-phase line with ASan/UBSan" OFF)
if(MPMC_SW92_THREE_PHASE_EXPERIMENTAL_LINE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 experimental-line sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_three_phase_experimental_line_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_three_phase_experimental_line_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

foreach(case IN ITEMS experimental_line permutation)
    add_test(NAME flash.sw92_three_phase_experimental_line.${case}
             COMMAND mpmc_sw92_three_phase_experimental_line_tests ${case})
    set_tests_properties(flash.sw92_three_phase_experimental_line.${case} PROPERTIES
        LABELS "flash;flash.sw92_three_phase_experimental_line;physical" TIMEOUT 180)
endforeach()
