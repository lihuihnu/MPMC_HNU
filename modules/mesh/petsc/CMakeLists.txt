cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_mesh_petsc LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_mesh)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/.." mesh-core)
endif()

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc=3.19.6)

add_library(mpmc_mesh_petsc INTERFACE)
add_library(mpmc::mesh_petsc ALIAS mpmc_mesh_petsc)
target_compile_features(mpmc_mesh_petsc INTERFACE cxx_std_20)
target_include_directories(
    mpmc_mesh_petsc
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_mesh_petsc
    INTERFACE mpmc::mesh MPI::MPI_CXX PkgConfig::PETSC)
