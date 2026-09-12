cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_split_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_CPA_SPLIT_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA split tests" OFF)

add_executable(mpmc_cpa_split_tests cpa_split_test.cpp)
target_link_libraries(mpmc_cpa_split_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_split_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_split_headers cpa_split_headers.cpp)
target_link_libraries(mpmc_cpa_split_headers PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_split_headers PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_exact_split_probe cpa_exact_split_probe.cpp)
target_link_libraries(mpmc_cpa_exact_split_probe PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_exact_split_probe PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    foreach(target IN ITEMS mpmc_cpa_split_tests mpmc_cpa_split_headers mpmc_cpa_exact_split_probe)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    endforeach()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(target IN ITEMS mpmc_cpa_split_tests mpmc_cpa_split_headers mpmc_cpa_exact_split_probe)
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endforeach()
endif()

if(MPMC_CPA_SPLIT_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA split sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN ITEMS mpmc_cpa_split_tests mpmc_cpa_split_headers mpmc_cpa_exact_split_probe)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

foreach(case IN ITEMS
        explicit_density_sides
        structural_anchor_residual
        associating_provider)
    add_test(NAME flash.cpa_split.${case}
             COMMAND mpmc_cpa_split_tests ${case})
    set_tests_properties(flash.cpa_split.${case} PROPERTIES
        LABELS "flash;cpa;split" TIMEOUT 300)
endforeach()

add_test(NAME flash.cpa_split.exact_probe COMMAND mpmc_cpa_exact_split_probe)
set_tests_properties(flash.cpa_split.exact_probe PROPERTIES
    LABELS "flash;cpa;split;diagnostic" TIMEOUT 300)

add_test(NAME flash.cpa_split.headers COMMAND mpmc_cpa_split_headers)
set_tests_properties(flash.cpa_split.headers PROPERTIES
    LABELS "flash;cpa;split;headers" TIMEOUT 60)
