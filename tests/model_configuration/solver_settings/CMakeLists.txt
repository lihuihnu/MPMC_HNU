cmake_minimum_required(VERSION 3.21)
project(mpmc_solver_settings_tests LANGUAGES CXX)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/model_configuration/pr76_solver" solver-library)
enable_testing()
add_executable(mpmc_solver_settings_tests settings_test.cpp solve_parity.cpp public_header.cpp)
target_link_libraries(mpmc_solver_settings_tests PRIVATE mpmc::pr76_solver_configuration)
set_target_properties(mpmc_solver_settings_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_solver_settings_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_solver_settings_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
option(MPMC_SOLVER_SETTINGS_ENABLE_SANITIZERS "Check immutable solver snapshots" OFF)
if(MPMC_SOLVER_SETTINGS_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "Solver settings sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_solver_settings_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_solver_settings_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()
foreach(case IN ITEMS preset custom_mapping required_fields invalid_domains resource_limits
        zero_false_semantics snapshot_isolation headers single_phase_parity three_phase_parity
        root_budget_parity initial_budget_parity final_two_budget_parity final_three_budget_parity)
    add_test(NAME model.solver_settings.${case} COMMAND mpmc_solver_settings_tests ${case})
    set_tests_properties(model.solver_settings.${case} PROPERTIES TIMEOUT 300 LABELS "model;solver_settings;parity")
endforeach()
