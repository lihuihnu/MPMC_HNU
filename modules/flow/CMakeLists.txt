cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_flow LANGUAGES CXX)
endif()

add_library(mpmc_flow INTERFACE)
add_library(mpmc::flow ALIAS mpmc_flow)
target_compile_features(mpmc_flow INTERFACE cxx_std_20)
target_include_directories(
    mpmc_flow
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
