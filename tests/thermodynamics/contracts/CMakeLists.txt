cmake_minimum_required(VERSION 3.21)
project(mpmc_thermodynamics_contract_tests LANGUAGES CXX)

# Independent consumer of the new module; no AD or root test project is configured.
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermo-library)
enable_testing()
add_executable(mpmc_thermodynamics_contract_tests
    contracts_test.cpp components_self_contained.cpp pr_self_contained.cpp)
target_link_libraries(mpmc_thermodynamics_contract_tests PRIVATE mpmc::thermodynamics)
set_target_properties(mpmc_thermodynamics_contract_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_thermodynamics_contract_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_thermodynamics_contract_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_THERMO_ENABLE_SANITIZERS "Check contract ownership with ASan/UBSan" OFF)
if(MPMC_THERMO_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "Contract sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_thermodynamics_contract_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_thermodynamics_contract_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(contract_cases ordered_identity order_mutations all_permutations identity_errors
    missing_parameters duplicate_and_pair_errors values_and_units provenance_policy model_revision
    applicability ownership_and_recovery resource_limits headers)
foreach(case IN LISTS contract_cases)
    add_test(NAME thermo.contracts.${case} COMMAND mpmc_thermodynamics_contract_tests ${case})
    set_tests_properties(thermo.contracts.${case} PROPERTIES LABELS "thermo;thermo.contracts;unit" TIMEOUT 30)
endforeach()
