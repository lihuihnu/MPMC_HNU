cmake_minimum_required(VERSION 3.21)
project(mpmc_ad_math_tests LANGUAGES CXX)

# This independent test project does not configure, compile, or run the old
# arithmetic suite. Only the AD INTERFACE target is consumed, as in a client.
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
enable_testing()

add_executable(mpmc_ad_math_tests math_test.cpp nested_math_test.cpp header_self_contained.cpp)
target_link_libraries(mpmc_ad_math_tests PRIVATE mpmc::ad)
set_target_properties(mpmc_ad_math_tests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(mpmc_ad_math_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_ad_math_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_AD_ENABLE_SANITIZERS "Check the new math tests with ASan and UBSan" OFF)
if(MPMC_AD_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "AD sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_ad_math_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_ad_math_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

set(math_cases
    exp exp2 expm1 log log2 log10 log1p sqrt cbrt
    sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh abs
    pow_integer pow_real pow_dual pow_scalar_base composition range_regressions
    overflow_policy nested_second_order header_odr)
foreach(case IN LISTS math_cases)
    add_test(NAME ad.math.${case} COMMAND mpmc_ad_math_tests ${case})
    set_tests_properties(ad.math.${case} PROPERTIES LABELS "ad;ad.math;unit" TIMEOUT 30)
endforeach()
