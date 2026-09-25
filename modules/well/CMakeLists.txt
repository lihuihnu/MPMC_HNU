cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_well LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_mesh)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../mesh"
        well-mesh-core)
endif()

add_library(mpmc_well INTERFACE)
add_library(mpmc::well ALIAS mpmc_well)
target_compile_features(mpmc_well INTERFACE cxx_std_20)
target_include_directories(
    mpmc_well
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_well
    INTERFACE mpmc::mesh)
