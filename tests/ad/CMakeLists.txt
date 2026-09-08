add_executable(mpmc_ad_tests dual_test.cpp header_self_contained.cpp)
target_link_libraries(mpmc_ad_tests PRIVATE mpmc::ad)
set_target_properties(mpmc_ad_tests PROPERTIES CXX_EXTENSIONS OFF)

# Diagnostics are private: a library consumer does not inherit test policies.
if(MSVC)
    target_compile_options(mpmc_ad_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_ad_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_AD_ENABLE_SANITIZERS "Use AddressSanitizer and UndefinedBehaviorSanitizer for AD tests" OFF)
if(MPMC_AD_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "The initial AD sanitizer configuration supports Linux GCC/Clang only")
    endif()
    target_compile_options(mpmc_ad_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_ad_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

add_test(NAME ad.dual COMMAND mpmc_ad_tests)
set_tests_properties(ad.dual PROPERTIES LABELS "ad;unit" TIMEOUT 30)
