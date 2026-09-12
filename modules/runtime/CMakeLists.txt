cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_runtime LANGUAGES CXX)
endif()

if(NOT TARGET mpmc::flash)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../flash" flash-library)
endif()

add_library(mpmc_runtime INTERFACE)
add_library(mpmc::runtime ALIAS mpmc_runtime)
target_compile_features(mpmc_runtime INTERFACE cxx_std_20)
target_include_directories(mpmc_runtime INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mpmc_runtime INTERFACE mpmc::flash)
