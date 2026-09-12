cmake_minimum_required(VERSION 3.21)
project(mpmc_pt_flash_backend_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_PT_FLASH_BACKEND_ENABLE_SANITIZERS
       "Enable ASan/UBSan for PT flash backend tests" OFF)

add_executable(mpmc_pt_flash_backend_tests
    pt_flash_backend_test.cpp
    pt_flash_backend_headers.cpp)
target_link_libraries(mpmc_pt_flash_backend_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_pt_flash_backend_tests PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../../thermodynamics/sw92")
set_target_properties(mpmc_pt_flash_backend_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_pt_flash_backend_tests PRIVATE
        /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_pt_flash_backend_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

if(MPMC_PT_FLASH_BACKEND_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PT flash backend sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_pt_flash_backend_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_pt_flash_backend_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(cases
    pr76_direct_equivalence
    sw92_direct_equivalence
    sw92_three_phase_capability
    contract_guards
    headers)
foreach(case IN LISTS cases)
    add_test(NAME flash.pt_flash_backend.${case}
             COMMAND mpmc_pt_flash_backend_tests ${case})
    set_tests_properties(flash.pt_flash_backend.${case} PROPERTIES
        LABELS "flash;pt_flash_backend;capability;conformance" TIMEOUT 300)
endforeach()
