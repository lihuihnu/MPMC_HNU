cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_physical_validation_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_CPA_PHYSICAL_VALIDATION_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA physical validation tests" OFF)

add_executable(mpmc_cpa_physical_validation_tests
    cpa_physical_validation_test.cpp)
target_link_libraries(mpmc_cpa_physical_validation_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_physical_validation_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_physical_validation_headers
    cpa_physical_validation_headers.cpp)
target_link_libraries(mpmc_cpa_physical_validation_headers PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_physical_validation_headers PROPERTIES CXX_EXTENSIONS OFF)

set(all_targets
    mpmc_cpa_physical_validation_tests
    mpmc_cpa_physical_validation_headers)

if(MSVC)
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    endforeach()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endforeach()
endif()

if(MPMC_CPA_PHYSICAL_VALIDATION_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA physical validation sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

foreach(case IN ITEMS
        parameter_provenance
        default_flash_precision_probe
        literature_vle_points
        component_permutation)
    add_test(NAME flash.cpa_physical_validation.${case}
             COMMAND mpmc_cpa_physical_validation_tests ${case})
    set_tests_properties(flash.cpa_physical_validation.${case} PROPERTIES
        LABELS "flash;cpa;physical_validation;association;literature" TIMEOUT 300)
endforeach()

add_test(NAME flash.cpa_physical_validation.headers
         COMMAND mpmc_cpa_physical_validation_headers)
set_tests_properties(flash.cpa_physical_validation.headers PROPERTIES
    LABELS "flash;cpa;physical_validation;headers" TIMEOUT 60)
