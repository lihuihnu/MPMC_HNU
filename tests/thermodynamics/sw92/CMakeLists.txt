cmake_minimum_required(VERSION 3.21)
project(mpmc_thermodynamics_sw92_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermo-library)
enable_testing()
add_executable(mpmc_thermodynamics_sw92_tests sw92_test.cpp headers.cpp)
target_link_libraries(mpmc_thermodynamics_sw92_tests PRIVATE mpmc::thermodynamics)
set_target_properties(mpmc_thermodynamics_sw92_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_thermodynamics_sw92_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_thermodynamics_sw92_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
option(MPMC_THERMO_ENABLE_SANITIZERS "Check SW92 integration with ASan/UBSan" OFF)
if(MPMC_THERMO_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_thermodynamics_sw92_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_thermodynamics_sw92_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(sw92_cases contract water_alpha bip_correlations mixing_reference phase_three_roots
    phase_aqueous_reference permutations background_pairs runtime_snapshots
    applicability_bounds input_domains headers)
foreach(case IN LISTS sw92_cases)
    add_test(NAME thermo.sw92.${case} COMMAND mpmc_thermodynamics_sw92_tests ${case})
    set_tests_properties(thermo.sw92.${case} PROPERTIES
        LABELS "thermo;thermo.sw92;unit" TIMEOUT 30)
endforeach()
