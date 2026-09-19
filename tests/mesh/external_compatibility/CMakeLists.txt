cmake_minimum_required(VERSION 3.21)
project(mpmc_mesh_external_compatibility LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/mesh"
    mesh-library)

add_executable(
    mpmc_mesh_external_compatibility
    external_mesh_compatibility.cpp)
target_link_libraries(
    mpmc_mesh_external_compatibility
    PRIVATE mpmc::mesh)
set_target_properties(
    mpmc_mesh_external_compatibility
    PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(
        mpmc_mesh_external_compatibility
        PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(
        mpmc_mesh_external_compatibility
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
            -Werror)
endif()
