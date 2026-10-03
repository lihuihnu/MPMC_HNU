cmake_minimum_required(VERSION 3.21)
project(mpmc_root_library_consumer LANGUAGES CXX)

# The parent owns its test; disabling repository tests must not disable this one.
enable_testing()
set(BUILD_TESTING OFF CACHE BOOL "Do not import the repository's AD tests" FORCE)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../.." mpmc-library)

if(NOT TARGET mpmc::ad OR TARGET mpmc_ad_tests)
    message(FATAL_ERROR "Root library-only build must expose AD without its tests")
endif()
foreach(module IN ITEMS thermodynamics flash mesh discretization
                        flow flow_thermodynamics flow_discretization well well_discretization)
    set(expected OFF)
    if(module STREQUAL "thermodynamics" AND
       (MPMC_ENABLE_THERMODYNAMICS OR MPMC_ENABLE_FLASH OR
        MPMC_ENABLE_FLOW_DISCRETIZATION OR MPMC_ENABLE_WELL_DISCRETIZATION))
        set(expected ON)
    elseif(module STREQUAL "flash" AND MPMC_ENABLE_FLASH)
        set(expected ON)
    elseif(module STREQUAL "mesh" AND
           (MPMC_ENABLE_MESH OR MPMC_ENABLE_DISCRETIZATION OR
            MPMC_ENABLE_FLOW_DISCRETIZATION OR MPMC_ENABLE_WELL_DISCRETIZATION))
        set(expected ON)
    elseif(module STREQUAL "discretization" AND
           (MPMC_ENABLE_DISCRETIZATION OR MPMC_ENABLE_FLOW_DISCRETIZATION OR MPMC_ENABLE_WELL_DISCRETIZATION))
        set(expected ON)
    elseif((module STREQUAL "flow" OR module STREQUAL "flow_thermodynamics" OR
            module STREQUAL "flow_discretization") AND
           (MPMC_ENABLE_FLOW_DISCRETIZATION OR MPMC_ENABLE_WELL_DISCRETIZATION))
        set(expected ON)
    elseif((module STREQUAL "well" OR module STREQUAL "well_discretization") AND
           MPMC_ENABLE_WELL_DISCRETIZATION)
        set(expected ON)
    endif()
    if(expected AND NOT TARGET mpmc::${module})
        message(FATAL_ERROR "Missing requested root library: ${module}")
    elseif(NOT expected AND TARGET mpmc::${module})
        message(FATAL_ERROR "Unrequested root library: ${module}")
    endif()
endforeach()

# Link only the highest requested target: all lower include paths and C++20
# requirements must arrive through production target usage requirements.
set(consumer_library ad)
if(MPMC_ENABLE_FLASH)
    set(consumer_library flash)
elseif(MPMC_ENABLE_THERMODYNAMICS)
    set(consumer_library thermodynamics)
endif()
add_executable(mpmc_root_library_consumer consumer.cpp)
target_link_libraries(mpmc_root_library_consumer PRIVATE mpmc::${consumer_library})
target_compile_definitions(mpmc_root_library_consumer PRIVATE
    MPMC_TEST_THERMODYNAMICS=$<OR:$<BOOL:${MPMC_ENABLE_THERMODYNAMICS}>,$<BOOL:${MPMC_ENABLE_FLASH}>>
    MPMC_TEST_FLASH=$<BOOL:${MPMC_ENABLE_FLASH}>)
set(consumer_targets mpmc_root_library_consumer)
if(MPMC_ENABLE_MESH OR MPMC_ENABLE_DISCRETIZATION)
    # A separate executable must obtain mesh only from the spatial chain,
    # without borrowing any include paths or dependencies from AD/thermo/flash.
    set(spatial_library mesh)
    if(MPMC_ENABLE_DISCRETIZATION)
        set(spatial_library discretization)
    endif()
    add_executable(mpmc_root_spatial_consumer spatial_consumer.cpp)
    target_link_libraries(mpmc_root_spatial_consumer PRIVATE mpmc::${spatial_library})
    target_compile_definitions(mpmc_root_spatial_consumer PRIVATE
        MPMC_TEST_DISCRETIZATION=$<BOOL:${MPMC_ENABLE_DISCRETIZATION}>)
    list(APPEND consumer_targets mpmc_root_spatial_consumer)
    add_test(NAME build.root_spatial_libraries COMMAND mpmc_root_spatial_consumer)
    set_tests_properties(build.root_spatial_libraries PROPERTIES LABELS "build;consumer" TIMEOUT 30)
endif()
if(MPMC_ENABLE_FLOW_DISCRETIZATION OR MPMC_ENABLE_WELL_DISCRETIZATION)
    # Keep the two published usage contracts separate: the discretization
    # bridge links flow + spatial targets, while flow_thermodynamics links EOS/AD.
    add_executable(mpmc_root_flow_consumer flow_consumer.cpp)
    target_link_libraries(mpmc_root_flow_consumer PRIVATE mpmc::flow_thermodynamics)
    add_executable(mpmc_root_flow_discretization_consumer flow_discretization_consumer.cpp)
    target_link_libraries(mpmc_root_flow_discretization_consumer PRIVATE mpmc::flow_discretization)
    list(APPEND consumer_targets mpmc_root_flow_consumer mpmc_root_flow_discretization_consumer)
    add_test(NAME build.root_flow_thermodynamics COMMAND mpmc_root_flow_consumer)
    add_test(NAME build.root_flow_discretization COMMAND mpmc_root_flow_discretization_consumer)
    set_tests_properties(build.root_flow_thermodynamics build.root_flow_discretization
                        PROPERTIES LABELS "build;consumer" TIMEOUT 30)
endif()
if(MPMC_ENABLE_WELL_DISCRETIZATION)
    # Check the base well target without borrowing the bridge's include paths.
    add_executable(mpmc_root_well_consumer well_consumer.cpp)
    target_link_libraries(mpmc_root_well_consumer PRIVATE mpmc::well)
    add_executable(mpmc_root_well_discretization_consumer well_discretization_consumer.cpp)
    target_link_libraries(mpmc_root_well_discretization_consumer PRIVATE mpmc::well_discretization)
    list(APPEND consumer_targets mpmc_root_well_consumer mpmc_root_well_discretization_consumer)
    add_test(NAME build.root_well COMMAND mpmc_root_well_consumer)
    add_test(NAME build.root_well_discretization COMMAND mpmc_root_well_discretization_consumer)
    set_tests_properties(build.root_well build.root_well_discretization
                        PROPERTIES LABELS "build;consumer" TIMEOUT 30)
endif()
if(MPMC_AD_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "Root consumer sanitizers require Linux GCC/Clang")
    endif()
endif()
foreach(consumer IN LISTS consumer_targets)
    set_target_properties(${consumer} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${consumer} PRIVATE /W4 /WX /permissive- /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(${consumer} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endif()
    if(MPMC_AD_ENABLE_SANITIZERS)
        target_compile_options(${consumer} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${consumer} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endif()
endforeach()
add_test(NAME build.root_libraries COMMAND mpmc_root_library_consumer)
set_tests_properties(build.root_libraries PROPERTIES LABELS "build;consumer" TIMEOUT 30)
