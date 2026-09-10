cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_flash LANGUAGES CXX)
endif()
if(NOT TARGET mpmc::thermodynamics)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../thermodynamics" thermo-library)
endif()
add_library(mpmc_flash INTERFACE)
add_library(mpmc::flash ALIAS mpmc_flash)
target_compile_features(mpmc_flash INTERFACE cxx_std_20)
target_include_directories(mpmc_flash INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mpmc_flash INTERFACE mpmc::thermodynamics)

# Optional converged-solution sensitivity adapter. Keep the established flash
# target free of an AD link dependency; consumers of pr76_sensitivity.hpp use
# mpmc::flash_sensitivity, which composes the independent AD module explicitly.
if(NOT TARGET mpmc::ad)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../ad" ad-library)
endif()
add_library(mpmc_flash_sensitivity INTERFACE)
add_library(mpmc::flash_sensitivity ALIAS mpmc_flash_sensitivity)
target_compile_features(mpmc_flash_sensitivity INTERFACE cxx_std_20)
target_link_libraries(mpmc_flash_sensitivity INTERFACE mpmc::flash mpmc::ad)
