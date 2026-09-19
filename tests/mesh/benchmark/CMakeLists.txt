cmake_minimum_required(VERSION 3.21)
project(mpmc_mesh_baseline_benchmark LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/mesh"
    mesh-library)

add_executable(
    mpmc_mesh_baseline_benchmark
    mesh_baseline_benchmark.cpp)
target_link_libraries(
    mpmc_mesh_baseline_benchmark
    PRIVATE mpmc::mesh)
set_target_properties(
    mpmc_mesh_baseline_benchmark
    PROPERTIES
        CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(
        mpmc_mesh_baseline_benchmark
        PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(
        mpmc_mesh_baseline_benchmark
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
            -Werror)
endif()
