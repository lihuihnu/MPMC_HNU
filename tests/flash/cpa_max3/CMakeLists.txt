cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_max3_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_CPA_MAX3_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA max3 tests" OFF)

add_executable(mpmc_cpa_max3_tests cpa_max3_test.cpp)
target_link_libraries(mpmc_cpa_max3_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_max3_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_max3_headers cpa_max3_headers.cpp)
target_link_libraries(mpmc_cpa_max3_headers PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_max3_headers PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    foreach(target IN ITEMS mpmc_cpa_max3_tests mpmc_cpa_max3_headers)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    endforeach()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(target IN ITEMS mpmc_cpa_max3_tests mpmc_cpa_max3_headers)
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endforeach()
endif()

if(MPMC_CPA_MAX3_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA max3 sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN ITEMS mpmc_cpa_max3_tests mpmc_cpa_max3_headers)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

foreach(case IN ITEMS
        fixed_three_phase
        two_to_three
        diagnostic_disappearance)
    add_test(NAME flash.cpa_max3.${case}
             COMMAND mpmc_cpa_max3_tests ${case})
    set_tests_properties(flash.cpa_max3.${case} PROPERTIES
        LABELS "flash;cpa;max3" TIMEOUT 180)
endforeach()

add_test(NAME flash.cpa_max3.headers COMMAND mpmc_cpa_max3_headers)
set_tests_properties(flash.cpa_max3.headers PROPERTIES
    LABELS "flash;cpa;max3;headers" TIMEOUT 60)
