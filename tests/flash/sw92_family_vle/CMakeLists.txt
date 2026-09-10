cmake_minimum_required(VERSION 3.21)
project(mpmc_sw92_family_vle_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

add_executable(mpmc_sw92_family_vle_tests
    sw92_family_vle_test.cpp
    sw92_split_header.cpp)
target_link_libraries(mpmc_sw92_family_vle_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_sw92_family_vle_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_sw92_family_vle_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_sw92_family_vle_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_SW92_FAMILY_VLE_ENABLE_SANITIZERS
       "Check SW92 one-family VLE with ASan/UBSan" OFF)
if(MPMC_SW92_FAMILY_VLE_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "SW92 family VLE sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_sw92_family_vle_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_sw92_family_vle_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    provider_roles
    co2_nonaqueous_vle
    methane_aqueous_vle
    permutation
    contracts_and_failures
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.sw92_family_vle.${case}
             COMMAND mpmc_sw92_family_vle_tests ${case})
    set_tests_properties(flash.sw92_family_vle.${case} PROPERTIES
        LABELS "flash;flash.sw92_family_vle;unit" TIMEOUT 120)
endforeach()
