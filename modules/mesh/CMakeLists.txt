cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_mesh LANGUAGES CXX)
endif()

add_library(mpmc_mesh INTERFACE)
add_library(mpmc::mesh ALIAS mpmc_mesh)
target_compile_features(mpmc_mesh INTERFACE cxx_std_20)
target_include_directories(mpmc_mesh INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
