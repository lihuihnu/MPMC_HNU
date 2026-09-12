cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_runtime_grpc LANGUAGES CXX)
endif()

if(NOT TARGET mpmc::runtime)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../runtime" runtime-library)
endif()

find_package(Protobuf REQUIRED)
find_package(gRPC CONFIG REQUIRED)

set(mpmc_pt_proto_root "${CMAKE_CURRENT_LIST_DIR}/../../api/proto")
set(mpmc_pt_proto
    "${mpmc_pt_proto_root}/mpmc/runtime/v1/pt_service.proto")
set(mpmc_pt_generated_root "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(mpmc_pt_generated_dir
    "${mpmc_pt_generated_root}/mpmc/runtime/v1")
set(mpmc_pt_generated_sources
    "${mpmc_pt_generated_dir}/pt_service.pb.cc"
    "${mpmc_pt_generated_dir}/pt_service.grpc.pb.cc")
set(mpmc_pt_generated_headers
    "${mpmc_pt_generated_dir}/pt_service.pb.h"
    "${mpmc_pt_generated_dir}/pt_service.grpc.pb.h")

if(TARGET protobuf::protoc)
    set(mpmc_protoc "$<TARGET_FILE:protobuf::protoc>")
elseif(Protobuf_PROTOC_EXECUTABLE)
    set(mpmc_protoc "${Protobuf_PROTOC_EXECUTABLE}")
else()
    message(FATAL_ERROR "A protoc executable is required for runtime_grpc")
endif()

if(TARGET gRPC::grpc_cpp_plugin)
    set(mpmc_grpc_cpp_plugin "$<TARGET_FILE:gRPC::grpc_cpp_plugin>")
else()
    find_program(mpmc_grpc_cpp_plugin grpc_cpp_plugin REQUIRED)
endif()

add_custom_command(
    OUTPUT ${mpmc_pt_generated_sources} ${mpmc_pt_generated_headers}
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${mpmc_pt_generated_dir}"
    COMMAND ${mpmc_protoc}
            "--proto_path=${mpmc_pt_proto_root}"
            "--cpp_out=${mpmc_pt_generated_root}"
            "--grpc_out=${mpmc_pt_generated_root}"
            "--plugin=protoc-gen-grpc=${mpmc_grpc_cpp_plugin}"
            "${mpmc_pt_proto}"
    DEPENDS "${mpmc_pt_proto}"
    COMMENT "Generating C++ Protobuf/gRPC sources for PT service v1"
    VERBATIM)

add_library(mpmc_pt_wire_generated STATIC
    ${mpmc_pt_generated_sources} ${mpmc_pt_generated_headers})
add_library(mpmc::pt_wire_generated ALIAS mpmc_pt_wire_generated)
target_compile_features(mpmc_pt_wire_generated PUBLIC cxx_std_20)
target_include_directories(mpmc_pt_wire_generated SYSTEM PUBLIC
    "${mpmc_pt_generated_root}")
target_link_libraries(mpmc_pt_wire_generated PUBLIC
    protobuf::libprotobuf gRPC::grpc++)
set_target_properties(mpmc_pt_wire_generated PROPERTIES CXX_EXTENSIONS OFF)

add_library(mpmc_runtime_grpc STATIC src/pt_grpc_adapter.cpp)
add_library(mpmc::runtime_grpc ALIAS mpmc_runtime_grpc)
add_dependencies(mpmc_runtime_grpc mpmc_pt_wire_generated)
target_compile_features(mpmc_runtime_grpc PUBLIC cxx_std_20)
target_include_directories(mpmc_runtime_grpc PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(mpmc_runtime_grpc PUBLIC
    mpmc::runtime mpmc::pt_wire_generated gRPC::grpc++)
set_target_properties(mpmc_runtime_grpc PROPERTIES CXX_EXTENSIONS OFF)
