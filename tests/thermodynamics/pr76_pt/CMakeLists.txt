cmake_minimum_required(VERSION 3.21)
project(mpmc_thermodynamics_pr76_pt_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/thermodynamics" thermo-library)
if(NOT TARGET mpmc_ad)
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
endif()
enable_testing()
add_executable(mpmc_thermodynamics_pr76_pt_tests pt_test.cpp phase_header.cpp roots_header.cpp)
add_executable(mpmc_thermodynamics_pr76_pt_ad_tests direct_ad_test.cpp)
target_link_libraries(mpmc_thermodynamics_pr76_pt_tests PRIVATE mpmc::thermodynamics mpmc::ad)
target_link_libraries(mpmc_thermodynamics_pr76_pt_ad_tests PRIVATE mpmc::thermodynamics mpmc::ad)
add_dependencies(mpmc_thermodynamics_pr76_pt_tests mpmc_thermodynamics_pr76_pt_ad_tests)
set_target_properties(mpmc_thermodynamics_pr76_pt_tests PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(mpmc_thermodynamics_pr76_pt_ad_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_thermodynamics_pr76_pt_tests PRIVATE /W4 /WX /permissive- /utf-8)
    target_compile_options(mpmc_thermodynamics_pr76_pt_ad_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_thermodynamics_pr76_pt_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    target_compile_options(mpmc_thermodynamics_pr76_pt_ad_tests PRIVATE
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
    target_compile_options(mpmc_thermodynamics_pr76_pt_ad_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_thermodynamics_pr76_pt_ad_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()
set(pr76_pt_cases roots_factored root_topology root_scale phase_values full_derivatives
    reduced_derivatives gibbs_duhem pure_limit low_pressure log_ratio zero_trace permutations
    runtime_shapes input_domains ownership_recovery attraction_cancellation derivative_guard headers)
foreach(case IN LISTS pr76_pt_cases)
    add_test(NAME thermo.pr76_pt.${case} COMMAND mpmc_thermodynamics_pr76_pt_tests ${case})
    set_tests_properties(thermo.pr76_pt.${case} PROPERTIES LABELS "thermo;thermo.pr76_pt;unit" TIMEOUT 30)
endforeach()
add_test(NAME thermo.pr76_pt.direct_ad COMMAND mpmc_thermodynamics_pr76_pt_ad_tests)
set_tests_properties(thermo.pr76_pt.direct_ad PROPERTIES LABELS "thermo;thermo.pr76_pt;unit" TIMEOUT 30)
