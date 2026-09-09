cmake_minimum_required(VERSION 3.21)
project(mpmc_flash_stability_tests LANGUAGES CXX)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
enable_testing()
add_executable(mpmc_flash_stability_tests stability_test.cpp generic_header.cpp pr_header.cpp)
target_link_libraries(mpmc_flash_stability_tests PRIVATE mpmc::flash mpmc::ad)
set_target_properties(mpmc_flash_stability_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_flash_stability_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_flash_stability_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()
option(MPMC_FLASH_ENABLE_SANITIZERS "Check stability integration with ASan/UBSan" OFF)
if(MPMC_FLASH_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "Stability sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_flash_stability_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_flash_stability_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()
set(cases ideal_distance ideal_search regular_distance regular_search trivial_is_not_proof
    zero_trace input_domains budgets failure_semantics backtracking_recovery roundoff_descent underflow_guard
    nonsmooth ownership_recovery
    minimum_gibbs_root pr_gradient pr_failure pr_permutations_shapes literature_nitrogen
    literature_carbon decimal_anchors headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.stability.${case} COMMAND mpmc_flash_stability_tests ${case})
    set_tests_properties(flash.stability.${case} PROPERTIES LABELS "flash;flash.stability;unit" TIMEOUT 60)
endforeach()
