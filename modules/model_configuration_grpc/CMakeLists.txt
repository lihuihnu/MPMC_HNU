cmake_minimum_required(VERSION 3.21)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(mpmc_model_configuration_grpc LANGUAGES CXX)
endif()
if(NOT TARGET mpmc::pr76_model_registry)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../model_configuration/pr76_registry" registry-library)
endif()
# Discover imported generator targets in this directory before adding children.
# Conan's protoc/plugin targets need not be GLOBAL or available on PATH.
find_package(Protobuf REQUIRED)
find_package(gRPC CONFIG REQUIRED)
# Reuse the frozen v1 wire data types, generation recipe and scoped gRPC header
# wrapper. The new service does not invoke or reinterpret the old PtService.
if(NOT TARGET mpmc::runtime_grpc)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../runtime_grpc" legacy-wire-library)
endif()
set(model_proto_root "${CMAKE_CURRENT_SOURCE_DIR}/proto")
set(legacy_proto_root "${CMAKE_CURRENT_SOURCE_DIR}/../../api/proto")
set(model_proto "${model_proto_root}/mpmc/model_configuration/v1/model_service.proto")
set(model_generated_root "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(model_generated_dir "${model_generated_root}/mpmc/model_configuration/v1")
set(model_sources "${model_generated_dir}/model_service.pb.cc" "${model_generated_dir}/model_service.grpc.pb.cc")
set(model_headers "${model_generated_dir}/model_service.pb.h" "${model_generated_dir}/model_service.grpc.pb.h")
if(TARGET protobuf::protoc)
    set(model_protoc "$<TARGET_FILE:protobuf::protoc>")
elseif(Protobuf_PROTOC_EXECUTABLE)
    set(model_protoc "${Protobuf_PROTOC_EXECUTABLE}")
else()
    find_program(model_protoc protoc REQUIRED)
endif()
if(TARGET gRPC::grpc_cpp_plugin)
    set(model_grpc_plugin "$<TARGET_FILE:gRPC::grpc_cpp_plugin>")
else()
    find_program(model_grpc_plugin grpc_cpp_plugin REQUIRED)
endif()
add_custom_command(OUTPUT ${model_sources} ${model_headers}
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${model_generated_dir}"
    COMMAND ${model_protoc} "--proto_path=${model_proto_root}" "--proto_path=${legacy_proto_root}"
        "--cpp_out=${model_generated_root}" "--grpc_out=${model_generated_root}"
        "--plugin=protoc-gen-grpc=${model_grpc_plugin}" "${model_proto}"
    DEPENDS "${model_proto}" "${legacy_proto_root}/mpmc/runtime/v1/pt_service.proto"
    VERBATIM)
add_library(mpmc_model_wire_generated STATIC ${model_sources} ${model_headers})
add_dependencies(mpmc_model_wire_generated mpmc_pt_wire_generated)
target_compile_features(mpmc_model_wire_generated PUBLIC cxx_std_20)
target_include_directories(mpmc_model_wire_generated SYSTEM PUBLIC "${model_generated_root}")
target_link_libraries(mpmc_model_wire_generated PUBLIC mpmc::pt_wire_generated)
if(MSVC)
    target_compile_options(mpmc_model_wire_generated PRIVATE
        "/FI${CMAKE_CURRENT_SOURCE_DIR}/../runtime_grpc/include/mpmc/runtime_grpc/grpc_headers.hpp")
endif()
add_library(mpmc_model_configuration_grpc STATIC src/model_wire.cpp src/model_grpc_adapter.cpp)
add_library(mpmc::model_configuration_grpc ALIAS mpmc_model_configuration_grpc)
add_dependencies(mpmc_model_configuration_grpc mpmc_model_wire_generated)
target_compile_features(mpmc_model_configuration_grpc PUBLIC cxx_std_20)
target_include_directories(mpmc_model_configuration_grpc PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../runtime_grpc/include")
target_link_libraries(mpmc_model_configuration_grpc PUBLIC mpmc::pr76_model_registry mpmc_model_wire_generated)
set_target_properties(mpmc_model_configuration_grpc mpmc_model_wire_generated PROPERTIES CXX_EXTENSIONS OFF)
