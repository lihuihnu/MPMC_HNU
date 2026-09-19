cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_discretization_petsc LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_mesh_petsc)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../mesh_petsc" mesh-petsc)
endif()

if(NOT TARGET mpmc_discretization)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../discretization"
        discretization-core)
endif()

add_library(mpmc_discretization_petsc INTERFACE)
add_library(mpmc::discretization_petsc ALIAS mpmc_discretization_petsc)
target_compile_features(mpmc_discretization_petsc INTERFACE cxx_std_20)
target_include_directories(
    mpmc_discretization_petsc
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_discretization_petsc
    INTERFACE mpmc::mesh_petsc mpmc::discretization)
