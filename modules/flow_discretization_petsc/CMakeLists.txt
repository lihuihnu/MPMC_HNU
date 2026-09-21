cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_flow_discretization_petsc LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_flow_discretization)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../flow_discretization"
        flow-discretization-core)
endif()

if(NOT TARGET mpmc_discretization_petsc)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../discretization_petsc"
        discretization-petsc)
endif()

if(NOT TARGET mpmc_flash)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../flash"
        flow-discretization-flash)
endif()

add_library(mpmc_flow_discretization_petsc INTERFACE)
add_library(
    mpmc::flow_discretization_petsc
    ALIAS mpmc_flow_discretization_petsc)
target_compile_features(
    mpmc_flow_discretization_petsc
    INTERFACE cxx_std_20)
target_include_directories(
    mpmc_flow_discretization_petsc
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_flow_discretization_petsc
    INTERFACE
        mpmc::flow_discretization
        mpmc::flow_thermodynamics
        mpmc::flash
        mpmc::discretization_petsc)
