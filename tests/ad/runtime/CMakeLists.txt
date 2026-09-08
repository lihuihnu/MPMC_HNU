cmake_minimum_required(VERSION 3.21)
project(mpmc_ad_runtime_tests LANGUAGES CXX)

# This independent test project does not configure, compile, or run the old
# arithmetic, elementary-function or fixed-Jacobian suites. Only the AD INTERFACE target
# is consumed, as in a client.
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
enable_testing()

add_executable(mpmc_ad_runtime_tests runtime_test.cpp header_self_contained.cpp)
target_link_libraries(mpmc_ad_runtime_tests PRIVATE mpmc::ad)
set_target_properties(mpmc_ad_runtime_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_ad_runtime_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_ad_runtime_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_AD_ENABLE_SANITIZERS "Check the new runtime tests with ASan and UBSan" OFF)
if(MPMC_AD_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "AD sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_ad_runtime_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_ad_runtime_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(runtime_cases
    block_seeding rectangular_analytic fixed_crosscheck constants_unused resize_ownership
    callback_contract failure_recovery input_validation size_limits output_validation
    primal_consistency workspace_isolation scaled_derivatives header_odr)
foreach(case IN LISTS runtime_cases)
    add_test(NAME ad.runtime.${case} COMMAND mpmc_ad_runtime_tests ${case})
    set_tests_properties(ad.runtime.${case} PROPERTIES LABELS "ad;ad.runtime;unit" TIMEOUT 30)
endforeach()
