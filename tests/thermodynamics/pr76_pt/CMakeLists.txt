cmake_minimum_required(VERSION 3.21)
project(mpmc_thermodynamics_pr76_pt_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermo-library)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
enable_testing()
add_executable(mpmc_thermodynamics_pr76_pt_tests pt_test.cpp phase_header.cpp roots_header.cpp)
target_link_libraries(mpmc_thermodynamics_pr76_pt_tests PRIVATE mpmc::thermodynamics mpmc::ad)
set_target_properties(mpmc_thermodynamics_pr76_pt_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_thermodynamics_pr76_pt_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_thermodynamics_pr76_pt_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
option(MPMC_THERMO_ENABLE_SANITIZERS "Check PT integration with ASan/UBSan" OFF)
if(MPMC_THERMO_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PR76 PT sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_thermodynamics_pr76_pt_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_thermodynamics_pr76_pt_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()
set(pr76_pt_cases roots_factored root_topology root_scale phase_values full_derivatives
    reduced_derivatives gibbs_duhem pure_limit low_pressure log_ratio zero_trace permutations
    runtime_shapes input_domains ownership_recovery attraction_cancellation derivative_guard headers)
foreach(case IN LISTS pr76_pt_cases)
    add_test(NAME thermo.pr76_pt.${case} COMMAND mpmc_thermodynamics_pr76_pt_tests ${case})
    set_tests_properties(thermo.pr76_pt.${case} PROPERTIES LABELS "thermo;thermo.pr76_pt;unit" TIMEOUT 30)
endforeach()
