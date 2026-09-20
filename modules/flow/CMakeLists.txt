cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_flow LANGUAGES CXX)
endif()

add_library(mpmc_flow INTERFACE)
add_library(mpmc::flow ALIAS mpmc_flow)
target_compile_features(mpmc_flow INTERFACE cxx_std_20)
target_include_directories(
    mpmc_flow
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")

if(NOT TARGET mpmc_thermodynamics)
    add_subdirectory(
        "${CMAKE_CURRENT_LIST_DIR}/../thermodynamics"
        flow-thermodynamics-library)
endif()

add_library(mpmc_flow_thermodynamics INTERFACE)
add_library(
    mpmc::flow_thermodynamics
    ALIAS mpmc_flow_thermodynamics)
target_compile_features(
    mpmc_flow_thermodynamics
    INTERFACE cxx_std_20)
target_link_libraries(
    mpmc_flow_thermodynamics
    INTERFACE
        mpmc::flow
        mpmc::thermodynamics)
