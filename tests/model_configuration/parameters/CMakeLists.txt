cmake_minimum_required(VERSION 3.21)
project(mpmc_model_parameter_tests LANGUAGES CXX)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/model_configuration" model-library)
enable_testing()
add_executable(mpmc_model_parameter_tests parameters_test.cpp public_header.cpp adapter_header.cpp)
target_link_libraries(mpmc_model_parameter_tests PRIVATE mpmc::model_configuration)
set_target_properties(mpmc_model_parameter_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_model_parameter_tests PRIVATE /W4 /WX /permissive- /utf-8)
else()
    target_compile_options(mpmc_model_parameter_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
option(MPMC_MODEL_ENABLE_SANITIZERS "Check parameter snapshot ownership" OFF)
if(MPMC_MODEL_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "Model sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_model_parameter_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_model_parameter_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()
foreach(case IN ITEMS complete_models order_mutations missing_parameters identity_pairs
        numeric_domains provenance applicability immutability quotas versions headers)
    add_test(NAME model.parameters.${case} COMMAND mpmc_model_parameter_tests ${case})
    set_tests_properties(model.parameters.${case} PROPERTIES TIMEOUT 30 LABELS "model;parameters;unit")
endforeach()
