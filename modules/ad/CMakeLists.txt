cmake_minimum_required(VERSION 3.21)
# This directory is also usable directly with add_subdirectory or cmake -S.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_ad LANGUAGES CXX)
endif()

add_library(mpmc_ad INTERFACE)
add_library(mpmc::ad ALIAS mpmc_ad)
target_compile_features(mpmc_ad INTERFACE cxx_std_20)
target_include_directories(mpmc_ad INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
# No project-module, third-party, test, network, or platform-specific dependencies.
