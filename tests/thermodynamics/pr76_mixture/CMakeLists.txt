cmake_minimum_required(VERSION 3.21)
project(mpmc_thermodynamics_pr76_mixture_tests LANGUAGES CXX)

# Independent new integration suite: no existing AD/contract test project is configured.
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermo-library)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
enable_testing()
add_executable(mpmc_thermodynamics_pr76_mixture_tests
    mixture_test.cpp header_self_contained.cpp)
target_link_libraries(mpmc_thermodynamics_pr76_mixture_tests PRIVATE mpmc::thermodynamics mpmc::ad)
set_target_properties(mpmc_thermodynamics_pr76_mixture_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_thermodynamics_pr76_mixture_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_thermodynamics_pr76_mixture_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_THERMO_ENABLE_SANITIZERS "Check new PR76 mixture integration with ASan/UBSan" OFF)
if(MPMC_THERMO_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PR76 mixture sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_thermodynamics_pr76_mixture_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_thermodynamics_pr76_mixture_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(pr76_mixture_cases printed_binary full_gradients reduced_gradients chain_rule
    zero_fractions permutations runtime_shapes normalization temperature_and_seeds pair_semantics
    principal_roots zero_attraction_float ownership_and_recovery numerical_scale headers)
foreach(case IN LISTS pr76_mixture_cases)
    add_test(NAME thermo.pr76_mixture.${case} COMMAND mpmc_thermodynamics_pr76_mixture_tests ${case})
    set_tests_properties(thermo.pr76_mixture.${case} PROPERTIES
        LABELS "thermo;thermo.pr76_mixture;unit" TIMEOUT 30)
endforeach()
