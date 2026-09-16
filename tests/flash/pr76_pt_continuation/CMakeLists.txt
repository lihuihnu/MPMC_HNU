cmake_minimum_required(VERSION 3.21)
project(mpmc_pr76_pt_continuation_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_PR76_PT_CONTINUATION_ENABLE_SANITIZERS
       "Enable ASan/UBSan for PR76 PT continuation tests" OFF)

add_executable(mpmc_pr76_pt_continuation_tests
    pr76_pt_continuation_test.cpp)
target_link_libraries(mpmc_pr76_pt_continuation_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_pt_continuation_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_pt_continuation_headers
    pr76_pt_continuation_headers.cpp)
target_link_libraries(mpmc_pr76_pt_continuation_headers PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_pt_continuation_headers PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_pr76_literature_reference_test
    pr76_literature_reference_test.cpp)
target_link_libraries(mpmc_pr76_literature_reference_test PRIVATE mpmc::flash)
set_target_properties(mpmc_pr76_literature_reference_test PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    foreach(target IN ITEMS
            mpmc_pr76_pt_continuation_tests
            mpmc_pr76_pt_continuation_headers
            mpmc_pr76_literature_reference_test)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    endforeach()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(target IN ITEMS
            mpmc_pr76_pt_continuation_tests
            mpmc_pr76_pt_continuation_headers
            mpmc_pr76_literature_reference_test)
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endforeach()
endif()

if(MPMC_PR76_PT_CONTINUATION_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "PR76 PT continuation sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN ITEMS
            mpmc_pr76_pt_continuation_tests
            mpmc_pr76_pt_continuation_headers
            mpmc_pr76_literature_reference_test)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

foreach(case IN ITEMS
        repeated_three_phase
        moving_three_phase
        invalid_path
        bracket_contract
        bidirectional_phase_sequence
        unresolved_resets)
    add_test(NAME flash.pr76_pt_continuation.${case}
             COMMAND mpmc_pr76_pt_continuation_tests ${case})
    set_tests_properties(flash.pr76_pt_continuation.${case} PROPERTIES
        LABELS "flash;pr76;continuation;max3" TIMEOUT 300)
endforeach()

add_test(NAME flash.pr76_pt_continuation.literature_reference
         COMMAND mpmc_pr76_literature_reference_test)
set_tests_properties(flash.pr76_pt_continuation.literature_reference PROPERTIES
    LABELS "flash;pr76;continuation;literature_reference" TIMEOUT 60)

add_test(NAME flash.pr76_pt_continuation.headers
         COMMAND mpmc_pr76_pt_continuation_headers)
set_tests_properties(flash.pr76_pt_continuation.headers PROPERTIES
    LABELS "flash;pr76;continuation;headers" TIMEOUT 60)
