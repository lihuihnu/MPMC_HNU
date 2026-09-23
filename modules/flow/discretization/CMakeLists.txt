cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_flow_discretization LANGUAGES CXX)
endif()

if(NOT TARGET mpmc_flow)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
        flow-core)
endif()

if(NOT TARGET mpmc_discretization)
    add_subdirectory(
        "${CMAKE_CURRENT_SOURCE_DIR}/../../discretization"
        discretization-core)
endif()

add_library(mpmc_flow_discretization INTERFACE)
add_library(
    mpmc::flow_discretization
    ALIAS mpmc_flow_discretization)
target_compile_features(
    mpmc_flow_discretization
    INTERFACE cxx_std_20)
target_include_directories(
    mpmc_flow_discretization
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
    mpmc_flow_discretization
    INTERFACE
        mpmc::flow
        mpmc::discretization)
