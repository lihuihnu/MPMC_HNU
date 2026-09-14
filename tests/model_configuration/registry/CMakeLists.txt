cmake_minimum_required(VERSION 3.21)
project(mpmc_registry_tests LANGUAGES CXX)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/model_configuration/pr76_registry" solver-library)
find_package(Threads REQUIRED)
enable_testing()
add_executable(mpmc_registry_tests registry_test.cpp public_header.cpp)
target_link_libraries(mpmc_registry_tests PRIVATE mpmc::pr76_model_registry Threads::Threads)
set_target_properties(mpmc_registry_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_registry_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_registry_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
option(MPMC_REGISTRY_ENABLE_SANITIZERS "Check registry ownership and capacity" OFF)
if(MPMC_REGISTRY_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "Registry sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_pr76_model_registry PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_compile_options(mpmc_registry_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_registry_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()
foreach(case IN ITEMS creation_parity handles registry_scope capacity configuration_failure
        parameter_quotas entropy_failure reserved_capacity close_during_creation
        release_during_solve registry_destruction close_lifetime lease_admission
        solve_exception parallel_models release_acquire_race concurrent_creation)
    add_test(NAME model.registry.${case} COMMAND mpmc_registry_tests ${case})
    set_tests_properties(model.registry.${case} PROPERTIES TIMEOUT 300 LABELS "model;registry;lifetime")
endforeach()
