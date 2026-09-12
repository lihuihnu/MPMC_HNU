cmake_minimum_required(VERSION 3.21)
project(mpmc_pt_process_host_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/pt_process"
                 pt-process-library)
enable_testing()

option(MPMC_PT_PROCESS_ENABLE_SANITIZERS
       "Enable ASan/UBSan for PT composition/host tests" OFF)
set(MPMC_PT_PROCESS_TEST_CERT_DIR "" CACHE PATH
    "Directory containing generated PT process mTLS test certificates")

add_executable(mpmc_pt_process_tests
    pt_process_test.cpp
    pt_process_headers.cpp)
target_link_libraries(mpmc_pt_process_tests PRIVATE mpmc::pt_process)
target_include_directories(mpmc_pt_process_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")

add_executable(mpmc_pt_secure_test_server secure_test_server.cpp)
target_link_libraries(mpmc_pt_secure_test_server PRIVATE mpmc::pt_process)

add_executable(mpmc_pt_secure_client secure_client.cpp)
target_link_libraries(mpmc_pt_secure_client PRIVATE mpmc::runtime_grpc)

set(all_targets
    mpmc_pt_process
    mpmc_pt_process_tests
    mpmc_pt_secure_test_server
    mpmc_pt_secure_client)
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

if(MPMC_PT_PROCESS_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PT process sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all
            -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

foreach(case IN ITEMS composition_root observer headers)
    add_test(NAME runtime.pt_process.${case}
             COMMAND mpmc_pt_process_tests ${case})
    set_tests_properties(runtime.pt_process.${case} PROPERTIES
        LABELS "runtime;grpc;process_host;composition" TIMEOUT 30)
endforeach()

if(MPMC_PT_PROCESS_TEST_CERT_DIR)
    add_test(NAME runtime.pt_process.mtls_host
        COMMAND mpmc_pt_process_tests mtls_host
            "${MPMC_PT_PROCESS_TEST_CERT_DIR}/backend.crt"
            "${MPMC_PT_PROCESS_TEST_CERT_DIR}/backend.key"
            "${MPMC_PT_PROCESS_TEST_CERT_DIR}/edge-ca.crt"
            "${MPMC_PT_PROCESS_TEST_CERT_DIR}/edge.crt"
            "${MPMC_PT_PROCESS_TEST_CERT_DIR}/edge.key"
            "${MPMC_PT_PROCESS_TEST_CERT_DIR}/backend-ca.crt"
            "${MPMC_PT_PROCESS_TEST_CERT_DIR}/browser.crt"
            "${MPMC_PT_PROCESS_TEST_CERT_DIR}/browser.key")
    set_tests_properties(runtime.pt_process.mtls_host PROPERTIES
        LABELS "runtime;grpc;process_host;tls;authentication" TIMEOUT 30)
endif()
