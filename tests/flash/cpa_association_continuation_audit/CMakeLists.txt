cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_association_continuation_audit LANGUAGES CXX)

set(MPMC_SOURCE_ROOT "" CACHE PATH "MPMC_HNU source root supplying production modules")
if(MPMC_SOURCE_ROOT STREQUAL "")
    get_filename_component(MPMC_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
endif()

add_subdirectory("${MPMC_SOURCE_ROOT}/modules/flash" flash-library)

add_executable(mpmc_cpa_association_continuation_contract contract_test.cpp)
target_link_libraries(mpmc_cpa_association_continuation_contract PRIVATE mpmc::flash)
target_include_directories(mpmc_cpa_association_continuation_contract PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../cpa_physical_validation")
set_target_properties(mpmc_cpa_association_continuation_contract PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_association_continuation_audit audit.cpp)
target_link_libraries(mpmc_cpa_association_continuation_audit PRIVATE mpmc::flash)
target_include_directories(mpmc_cpa_association_continuation_audit PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../cpa_physical_validation")
set_target_properties(mpmc_cpa_association_continuation_audit PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_association_certificate_audit certificate_audit.cpp)
target_link_libraries(mpmc_cpa_association_certificate_audit PRIVATE mpmc::flash)
target_include_directories(mpmc_cpa_association_certificate_audit PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/../cpa_physical_validation")
set_target_properties(mpmc_cpa_association_certificate_audit PROPERTIES CXX_EXTENSIONS OFF)

foreach(target IN ITEMS
        mpmc_cpa_association_continuation_contract
        mpmc_cpa_association_continuation_audit
        mpmc_cpa_association_certificate_audit)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endif()
endforeach()
