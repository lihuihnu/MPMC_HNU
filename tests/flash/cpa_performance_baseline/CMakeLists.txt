cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_performance_baseline LANGUAGES CXX)

set(MPMC_SOURCE_ROOT "" CACHE PATH "MPMC_HNU source root supplying production modules")
if(MPMC_SOURCE_ROOT STREQUAL "")
    get_filename_component(MPMC_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
endif()

option(MPMC_CPA_PERF_ENABLE_COVERAGE
       "Enable GCC gcov instrumentation for structural hot-path counts" OFF)

add_subdirectory("${MPMC_SOURCE_ROOT}/modules/flash" flash-library)

add_executable(mpmc_cpa_performance_baseline benchmark.cpp)
target_link_libraries(mpmc_cpa_performance_baseline PRIVATE mpmc::flash)
target_include_directories(mpmc_cpa_performance_baseline PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../cpa_physical_validation")
set_target_properties(mpmc_cpa_performance_baseline PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_cpa_performance_baseline PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_cpa_performance_baseline PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

if(MPMC_CPA_PERF_ENABLE_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "CPA performance coverage counting requires GCC/gcov")
    endif()
    target_compile_options(mpmc_cpa_performance_baseline PRIVATE
        --coverage -fno-omit-frame-pointer)
    target_link_options(mpmc_cpa_performance_baseline PRIVATE --coverage)
endif()
