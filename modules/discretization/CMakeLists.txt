cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_discretization LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_mesh)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../mesh" mesh-core)
endif()

add_library(mpmc_discretization INTERFACE)
add_library(mpmc::discretization ALIAS mpmc_discretization)
target_compile_features(mpmc_discretization INTERFACE cxx_std_20)
target_include_directories(
    mpmc_discretization
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_discretization
    INTERFACE mpmc::mesh)
