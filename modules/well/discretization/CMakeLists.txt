cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_well_discretization LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_well)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
        well-core)
endif()

if(NOT TARGET mpmc_flow)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../../flow"
        well-discretization-flow-core)
endif()

if(NOT TARGET mpmc_flow_discretization)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../../flow/discretization"
        well-discretization-flow-discretization)
endif()

add_library(mpmc_well_discretization INTERFACE)
add_library(
    mpmc::well_discretization
    ALIAS mpmc_well_discretization)
target_compile_features(
    mpmc_well_discretization
    INTERFACE cxx_std_20)
target_include_directories(
    mpmc_well_discretization
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_well_discretization
    INTERFACE
        mpmc::well
        mpmc::flow
        mpmc::flow_discretization)
