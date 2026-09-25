cmake_minimum_required(VERSION 3.21)
project(mpmc_ad_jacobian_tests LANGUAGES CXX)

# This independent test project does not configure, compile, or run the old
# arithmetic or elementary-function suites. Only the AD INTERFACE target
# is consumed, as in a client.
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
enable_testing()

add_executable(mpmc_ad_jacobian_tests jacobian_test.cpp nested_hessian_test.cpp header_self_contained.cpp)
target_link_libraries(mpmc_ad_jacobian_tests PRIVATE mpmc::ad)
set_target_properties(mpmc_ad_jacobian_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_ad_jacobian_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_ad_jacobian_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_AD_ENABLE_SANITIZERS "Check the new jacobian tests with ASan and UBSan" OFF)
if(MPMC_AD_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "AD sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_ad_jacobian_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_ad_jacobian_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(jacobian_cases
    identity_seeds rectangular_analytic constants_unused_inputs single_input_output
    callable_forwarding repeated_calls_ownership exception_propagation
    nonfinite_inputs output_passthrough nested_hessian header_odr)
foreach(case IN LISTS jacobian_cases)
    add_test(NAME ad.jacobian.${case} COMMAND mpmc_ad_jacobian_tests ${case})
    set_tests_properties(ad.jacobian.${case} PROPERTIES LABELS "ad;ad.jacobian;unit" TIMEOUT 30)
endforeach()
