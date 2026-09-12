cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_pt_process LANGUAGES CXX)
endif()

if(NOT TARGET mpmc::runtime_grpc)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../runtime_grpc"
                     runtime-grpc-library)
endif()

add_library(mpmc_pt_process STATIC
    src/composition_root.cpp
    src/configured_backends.cpp
    src/json_line_observer.cpp
    src/process_host.cpp)
add_library(mpmc::pt_process ALIAS mpmc_pt_process)
target_compile_features(mpmc_pt_process PUBLIC cxx_std_20)
target_include_directories(mpmc_pt_process PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mpmc_pt_process PUBLIC mpmc::runtime_grpc)
set_target_properties(mpmc_pt_process PROPERTIES CXX_EXTENSIONS OFF)
