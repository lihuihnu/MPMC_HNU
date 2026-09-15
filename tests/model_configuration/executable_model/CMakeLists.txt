cmake_minimum_required(VERSION 3.21)
project(mpmc_executable_model_tests LANGUAGES CXX)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/model_configuration/pr76_solver" solver-library)
find_package(Threads REQUIRED)
enable_testing()
add_executable(mpmc_executable_model_tests model_test.cpp public_header.cpp)
add_executable(mpmc_executable_applicability_tests applicability_endpoint_test.cpp)
add_executable(mpmc_executable_literature_three_phase_test literature_three_phase_configured_test.cpp)
target_link_libraries(mpmc_executable_model_tests PRIVATE mpmc::pr76_solver_configuration Threads::Threads)
target_link_libraries(mpmc_executable_applicability_tests PRIVATE mpmc::pr76_solver_configuration Threads::Threads)
target_link_libraries(mpmc_executable_literature_three_phase_test PRIVATE mpmc::pr76_solver_configuration Threads::Threads)
option(MPMC_EXECUTABLE_MODEL_ENABLE_SANITIZERS "Check executable model ownership" OFF)
foreach(target IN ITEMS mpmc_executable_model_tests mpmc_executable_applicability_tests
        mpmc_executable_literature_three_phase_test)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endif()
    if(MPMC_EXECUTABLE_MODEL_ENABLE_SANITIZERS)
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
            message(FATAL_ERROR "Executable model sanitizers require Linux GCC/Clang")
        endif()
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endif()
endforeach()
foreach(case IN ITEMS single_phase two_phase ternary_cold_start parameter_isolation
        root_budget settings_isolation lifetime creation_errors resource_limits
        request_errors request_locations applicability public_hint_parity
        continuation_hint_parity hint_validation admission parallel_models)
    add_test(NAME model.executable.${case} COMMAND mpmc_executable_model_tests ${case})
    set_tests_properties(model.executable.${case} PROPERTIES TIMEOUT 300 LABELS "model;executable;parity")
endforeach()
add_test(NAME model.executable.applicability_endpoints COMMAND mpmc_executable_applicability_tests)
set_tests_properties(model.executable.applicability_endpoints PROPERTIES
    TIMEOUT 300 LABELS "model;executable;parity;applicability")
add_test(NAME model.executable.literature_three_phase_configured
         COMMAND mpmc_executable_literature_three_phase_test)
set_tests_properties(model.executable.literature_three_phase_configured PROPERTIES
    TIMEOUT 300 LABELS "model;executable;parity;literature;three_phase")
