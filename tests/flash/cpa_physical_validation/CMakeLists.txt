cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_physical_validation_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/flash" flash-library)
enable_testing()

option(MPMC_CPA_PHYSICAL_VALIDATION_ENABLE_SANITIZERS
       "Enable ASan/UBSan for CPA physical validation tests" OFF)

set(thermopack_oracle_file
    "${CMAKE_CURRENT_LIST_DIR}/../cpa_thermopack_oracle/thermopack_d68c794_meoh_h2o_33315k.json")
file(READ "${thermopack_oracle_file}" thermopack_oracle_json)
string(JSON thermopack_oracle_schema GET "${thermopack_oracle_json}" schema)
string(JSON thermopack_oracle_commit GET "${thermopack_oracle_json}" generator commit)
string(JSON thermopack_oracle_eos GET "${thermopack_oracle_json}" model physical_eos)
string(JSON thermopack_oracle_association_rule
       GET "${thermopack_oracle_json}" model association_delta_rule)
string(JSON thermopack_oracle_state_count LENGTH "${thermopack_oracle_json}" states)
if(NOT thermopack_oracle_schema STREQUAL
       "MPMC_HNU/CPA/ThermoPack-two-phase-TP-oracle/v1")
    message(FATAL_ERROR "unexpected frozen ThermoPack CPA oracle schema")
endif()
if(NOT thermopack_oracle_commit STREQUAL
       "d68c794c7342bfc6938eb424a1fbb88b7780b738")
    message(FATAL_ERROR "unexpected frozen ThermoPack revision")
endif()
if(NOT thermopack_oracle_eos STREQUAL "SRK" OR
   NOT thermopack_oracle_association_rule STREQUAL "STANDARD/CR-1")
    message(FATAL_ERROR "frozen ThermoPack oracle no longer matches the audited CPA formulation")
endif()
if(NOT thermopack_oracle_state_count EQUAL 5)
    message(FATAL_ERROR "frozen ThermoPack CPA oracle must contain exactly five states")
endif()

set(thermopack_state_initializers "")
foreach(index RANGE 0 4)
    string(JSON state_temperature GET "${thermopack_oracle_json}" states ${index} temperature_k)
    string(JSON state_pressure GET "${thermopack_oracle_json}" states ${index} pressure_pa)
    string(JSON state_feed_meoh GET "${thermopack_oracle_json}" states ${index} feed MEOH)
    string(JSON state_phase_code GET "${thermopack_oracle_json}" states ${index} phase_code)
    string(JSON state_beta_vapor GET "${thermopack_oracle_json}" states ${index} beta_vapor)
    string(JSON state_beta_liquid GET "${thermopack_oracle_json}" states ${index} beta_liquid)
    string(JSON state_liquid_meoh GET "${thermopack_oracle_json}" states ${index} liquid MEOH)
    string(JSON state_vapor_meoh GET "${thermopack_oracle_json}" states ${index} vapor MEOH)
    string(APPEND thermopack_state_initializers
        "    ThermoPackState{${state_temperature}, ${state_pressure}, ${state_feed_meoh}, ${state_phase_code}, ${state_beta_vapor}, ${state_beta_liquid}, ${state_liquid_meoh}, ${state_vapor_meoh}},\n")
endforeach()

set(thermopack_generated_header
    "${CMAKE_CURRENT_BINARY_DIR}/cpa_thermopack_oracle_generated.hpp")
file(GENERATE OUTPUT "${thermopack_generated_header}" CONTENT
"#ifndef MPMC_TEST_CPA_THERMOPACK_ORACLE_GENERATED_HPP\n#define MPMC_TEST_CPA_THERMOPACK_ORACLE_GENERATED_HPP\n\n#include <array>\n\nnamespace cpa_thermopack_oracle {\nstruct ThermoPackState {\n    double temperature_k;\n    double pressure_pa;\n    double feed_methanol;\n    int phase_code;\n    double beta_vapor;\n    double beta_liquid;\n    double liquid_methanol;\n    double vapor_methanol;\n};\n\ninline constexpr std::array<ThermoPackState, ${thermopack_oracle_state_count}> states{{\n${thermopack_state_initializers}}};\n\n} // namespace cpa_thermopack_oracle\n\n#endif // MPMC_TEST_CPA_THERMOPACK_ORACLE_GENERATED_HPP\n")

add_executable(mpmc_cpa_physical_validation_tests
    cpa_physical_validation_test.cpp)
target_link_libraries(mpmc_cpa_physical_validation_tests PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_physical_validation_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_physical_validation_headers
    cpa_physical_validation_headers.cpp)
target_link_libraries(mpmc_cpa_physical_validation_headers PRIVATE mpmc::flash)
set_target_properties(mpmc_cpa_physical_validation_headers PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_thermopack_delta_tests
    cpa_thermopack_delta_test.cpp)
target_link_libraries(mpmc_cpa_thermopack_delta_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_cpa_thermopack_delta_tests PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}")
set_target_properties(mpmc_cpa_thermopack_delta_tests PROPERTIES CXX_EXTENSIONS OFF)

set(all_targets
    mpmc_cpa_physical_validation_tests
    mpmc_cpa_physical_validation_headers
    mpmc_cpa_thermopack_delta_tests)

if(MSVC)
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
    endforeach()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endforeach()
endif()

if(MPMC_CPA_PHYSICAL_VALIDATION_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "CPA physical validation sanitizers require Linux GCC/Clang")
    endif()
    foreach(target IN LISTS all_targets)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-sanitize-recover=all)
    endforeach()
endif()

foreach(case IN ITEMS
        parameter_provenance
        default_flash_precision
        final_tpd_56652_safety)
    add_test(NAME flash.cpa_physical_validation.${case}
             COMMAND mpmc_cpa_physical_validation_tests ${case})
    set_tests_properties(flash.cpa_physical_validation.${case} PROPERTIES
        LABELS "flash;cpa;physical_validation;association;literature;full_physical;sanitizer_safety"
        TIMEOUT 300)
endforeach()

foreach(case IN ITEMS
        literature_vle_points
        component_permutation)
    add_test(NAME flash.cpa_physical_validation.${case}
             COMMAND mpmc_cpa_physical_validation_tests ${case})
    set_tests_properties(flash.cpa_physical_validation.${case} PROPERTIES
        LABELS "flash;cpa;physical_validation;association;literature;full_physical"
        TIMEOUT 300)
endforeach()

add_test(NAME flash.cpa_physical_validation.thermopack_five_state_delta
         COMMAND mpmc_cpa_thermopack_delta_tests)
set_tests_properties(flash.cpa_physical_validation.thermopack_five_state_delta PROPERTIES
    LABELS "flash;cpa;physical_validation;thermopack;external_oracle;full_physical"
    TIMEOUT 300)

add_test(NAME flash.cpa_physical_validation.headers
         COMMAND mpmc_cpa_physical_validation_headers)
set_tests_properties(flash.cpa_physical_validation.headers PROPERTIES
    LABELS "flash;cpa;physical_validation;headers;full_physical;sanitizer_safety"
    TIMEOUT 60)
