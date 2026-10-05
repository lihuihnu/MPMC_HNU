cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_simulation LANGUAGES CXX)
endif()

add_library(mpmc_simulation INTERFACE)
add_library(mpmc::simulation ALIAS mpmc_simulation)
target_compile_features(mpmc_simulation INTERFACE cxx_std_20)
target_include_directories(
    mpmc_simulation
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
# The lifecycle core is solver- and platform-independent: standard library only.
