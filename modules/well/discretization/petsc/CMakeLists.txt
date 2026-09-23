cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_well_discretization_petsc LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_well_discretization)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
        well-discretization-core)
endif()

if(NOT TARGET mpmc_flow_discretization_petsc)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../flow/discretization/petsc"
        well-discretization-flow-petsc)
endif()

add_library(mpmc_well_discretization_petsc INTERFACE)
add_library(
    mpmc::well_discretization_petsc
    ALIAS mpmc_well_discretization_petsc)
target_compile_features(
    mpmc_well_discretization_petsc
    INTERFACE cxx_std_20)
target_include_directories(
    mpmc_well_discretization_petsc
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_well_discretization_petsc
    INTERFACE
        mpmc::well_discretization
        mpmc::flow_discretization_petsc)
