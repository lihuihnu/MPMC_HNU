cmake_minimum_required(VERSION 3.21)
project(mpmc_pt_grpc_adapter_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/runtime_grpc"
                 runtime-grpc-library)
enable_testing()

option(MPMC_PT_GRPC_ENABLE_SANITIZERS
       "Enable ASan/UBSan for PT gRPC adapter tests" OFF)

add_executable(mpmc_pt_grpc_adapter_tests
    pt_grpc_adapter_test.cpp
    pt_grpc_adapter_headers.cpp)
target_link_libraries(mpmc_pt_grpc_adapter_tests PRIVATE mpmc::runtime_grpc)

add_executable(mpmc_pt_grpc_golden_server golden_server.cpp)
target_link_libraries(mpmc_pt_grpc_golden_server PRIVATE mpmc::runtime_grpc)

set(all_targets
    mpmc_runtime_grpc
    mpmc_pt_grpc_adapter_tests
    mpmc_pt_grpc_golden_server)
foreach(target IN LISTS all_targets)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
endforeach()

if(MSVC)
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    endforeach()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endforeach()
endif()

if(MPMC_PT_GRPC_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PT gRPC sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all
            -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

set(cases
    mapping_and_outcomes
    wire_and_size_limits
    bearer_authentication
    concurrency_gate
    deadline_and_cancel
    headers)
foreach(case IN LISTS cases)
    add_test(NAME runtime.pt_grpc_adapter.${case}
             COMMAND mpmc_pt_grpc_adapter_tests ${case})
    set_tests_properties(runtime.pt_grpc_adapter.${case} PROPERTIES
        LABELS "runtime;grpc;process_adapter;contract" TIMEOUT 20)
endforeach()
