cmake_minimum_required(VERSION 3.21)
project(mpmc_well_core_tests LANGUAGES CXX)

if(NOT TARGET mpmc_well)
    add_subdirectory(
        "${CMAKE_CURRENT_LIST_DIR}/../../../modules/well"
        well-library)
endif()
enable_testing()

option(
    MPMC_WELL_ENABLE_SANITIZERS
    "Use AddressSanitizer and UndefinedBehaviorSanitizer for well core tests"
    OFF)

function(configure_well_test_target target)
    target_link_libraries(
        ${target}
        PRIVATE mpmc::well)
    set_target_properties(
        ${target}
        PROPERTIES CXX_EXTENSIONS OFF)

    if(MSVC)
        target_compile_options(
            ${target}
            PRIVATE /W4 /WX /permissive- /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(
            ${target}
            PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endif()

    if(MPMC_WELL_ENABLE_SANITIZERS)
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
           NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
            message(FATAL_ERROR
                "Well core sanitizers require Linux GCC/Clang")
        endif()
        target_compile_options(
            ${target}
            PRIVATE
                -fsanitize=address,undefined
                -fno-sanitize-recover=all
                -fno-omit-frame-pointer)
        target_link_options(
            ${target}
            PRIVATE
                -fsanitize=address,undefined
                -fno-sanitize-recover=all)
    endif()
endfunction()

add_executable(
    mpmc_well_peaceman_tests
    peaceman_well_index_3d_test.cpp)
configure_well_test_target(
    mpmc_well_peaceman_tests)

add_executable(
    mpmc_well_peaceman_header_self_contained
    peaceman_well_index_3d_header_self_contained.cpp)
configure_well_test_target(
    mpmc_well_peaceman_header_self_contained)

add_custom_target(
    mpmc_well_core_all
    DEPENDS
        mpmc_well_peaceman_tests
        mpmc_well_peaceman_header_self_contained)

foreach(case IN ITEMS
    isotropic_square_limit
    anisotropic_rectangular_z
    axis_permutation_and_cell_binding
    skin_factor
    invalid_inputs)
    add_test(
        NAME discretization.core.well.peaceman_well_index_3d.${case}
        COMMAND mpmc_well_peaceman_tests ${case})
    set_tests_properties(
        discretization.core.well.peaceman_well_index_3d.${case}
        PROPERTIES
            LABELS "well;well.core;peaceman;well_index;unit"
            TIMEOUT 30)
endforeach()

add_test(
    NAME discretization.core.well.peaceman_well_index_3d.header_self_contained
    COMMAND mpmc_well_peaceman_header_self_contained)
set_tests_properties(
    discretization.core.well.peaceman_well_index_3d.header_self_contained
    PROPERTIES
        LABELS "well;well.core;peaceman;header"
        TIMEOUT 30)
