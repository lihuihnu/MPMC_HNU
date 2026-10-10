cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_simulation_petsc LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_simulation)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
        simulation-core)
endif()

if(NOT TARGET mpmc_flow_discretization_petsc)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../../flow/discretization/petsc"
        simulation-flow-discretization-petsc)
endif()

add_library(mpmc_simulation_petsc INTERFACE)
add_library(mpmc::simulation_petsc ALIAS mpmc_simulation_petsc)
target_compile_features(mpmc_simulation_petsc INTERFACE cxx_std_20)
target_include_directories(
    mpmc_simulation_petsc
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_simulation_petsc
    INTERFACE
        mpmc::simulation
        mpmc::flow_discretization_petsc)
