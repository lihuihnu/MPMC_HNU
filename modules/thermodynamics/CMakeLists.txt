cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_thermodynamics LANGUAGES CXX)
endif()

if(NOT TARGET mpmc::ad)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../ad" ad-library)
endif()

add_library(mpmc_thermodynamics INTERFACE)
add_library(mpmc::thermodynamics ALIAS mpmc_thermodynamics)
target_compile_features(mpmc_thermodynamics INTERFACE cxx_std_20)
target_include_directories(mpmc_thermodynamics INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mpmc_thermodynamics INTERFACE mpmc::ad)
