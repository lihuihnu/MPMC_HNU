cmake_minimum_required(VERSION 3.21)
project(mpmc_pr76_three_phase_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_PR76_THREE_PHASE_ENABLE_SANITIZERS
       "Enable ASan/UBSan for PR76 three-phase tests" OFF)

add_executable(mpmc_pr76_three_phase_tests
    pr76_three_phase_test.cpp
    pr76_three_phase_headers.cpp)
target_link_libraries(mpmc_pr76_three_phase_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_three_phase_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_pr76_three_phase_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_pr76_three_phase_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

if(MPMC_PR76_THREE_PHASE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PR76 three-phase sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_pr76_three_phase_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_pr76_three_phase_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    generic_exact
    generic_disappearance
    pr76_synthetic_max3
    pr76_fresh_two_phase_neighbor
    pr76_component_permutation
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.pr76_three_phase.${case}
             COMMAND mpmc_pr76_three_phase_tests ${case})
    set_tests_properties(flash.pr76_three_phase.${case} PROPERTIES
        LABELS "flash;pr76;three_phase;max3" TIMEOUT 300)
endforeach()
