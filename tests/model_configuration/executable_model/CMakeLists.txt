cmake_minimum_required(VERSION 3.21)
project(mpmc_executable_model_tests LANGUAGES CXX)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/model_configuration/pr76_solver" solver-library)
find_package(Threads REQUIRED)
enable_testing()
add_executable(mpmc_executable_model_tests model_test.cpp public_header.cpp)
target_link_libraries(mpmc_executable_model_tests PRIVATE mpmc::pr76_solver_configuration Threads::Threads)
set_target_properties(mpmc_executable_model_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_executable_model_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_executable_model_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
option(MPMC_EXECUTABLE_MODEL_ENABLE_SANITIZERS "Check executable model ownership" OFF)
if(MPMC_EXECUTABLE_MODEL_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "Executable model sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_executable_model_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_executable_model_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()
foreach(case IN ITEMS single_phase two_phase ternary_cold_start parameter_isolation
        root_budget settings_isolation lifetime creation_errors resource_limits
        request_errors request_locations applicability admission parallel_models)
    add_test(NAME model.executable.${case} COMMAND mpmc_executable_model_tests ${case})
    set_tests_properties(model.executable.${case} PROPERTIES TIMEOUT 300 LABELS "model;executable;parity")
endforeach()
