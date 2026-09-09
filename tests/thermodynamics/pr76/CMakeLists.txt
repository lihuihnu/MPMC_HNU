cmake_minimum_required(VERSION 3.21)
project(mpmc_thermodynamics_pr76_tests LANGUAGES CXX)

# Independent new integration suite: no existing AD/contract test project is configured.
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermo-library)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
enable_testing()
add_executable(mpmc_thermodynamics_pr76_tests
    pr76_test.cpp header_self_contained.cpp)
target_link_libraries(mpmc_thermodynamics_pr76_tests PRIVATE mpmc::thermodynamics mpmc::ad)
set_target_properties(mpmc_thermodynamics_pr76_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_thermodynamics_pr76_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_thermodynamics_pr76_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_THERMO_ENABLE_SANITIZERS "Check new PR76 integration with ASan/UBSan" OFF)
if(MPMC_THERMO_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PR76 sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_thermodynamics_pr76_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_thermodynamics_pr76_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(pr76_cases printed_reference analytic_temperature chain_rule fixed_jacobian runtime_jacobian
    snapshot_lifetime temperature_domain declared_bounds algebraic_extension numerical_range headers)
foreach(case IN LISTS pr76_cases)
    add_test(NAME thermo.pr76.${case} COMMAND mpmc_thermodynamics_pr76_tests ${case})
    set_tests_properties(thermo.pr76.${case} PROPERTIES LABELS "thermo;thermo.pr76;unit" TIMEOUT 30)
endforeach()
