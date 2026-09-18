cmake_minimum_required(VERSION 3.21)
project(mpmc_cpa_physical_validation_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/ad" ad-library)
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

set(thermopack_phase_kernel_file
    "${CMAKE_CURRENT_LIST_DIR}/../cpa_thermopack_oracle/thermopack_d68c794_meoh_h2o_phase_kernel_33315k.json")
file(READ "${thermopack_phase_kernel_file}" thermopack_phase_kernel_json)
string(JSON thermopack_phase_kernel_schema GET "${thermopack_phase_kernel_json}" schema)
string(JSON thermopack_phase_kernel_commit GET "${thermopack_phase_kernel_json}" source commit)
string(JSON thermopack_phase_kernel_rgas GET "${thermopack_phase_kernel_json}" gas_constant_j_per_mol_k)
string(JSON thermopack_phase_kernel_tc_meoh
       GET "${thermopack_phase_kernel_json}" alpha_critical_temperature_source MEOH temperature_k)
string(JSON thermopack_phase_kernel_tc_h2o
       GET "${thermopack_phase_kernel_json}" alpha_critical_temperature_source H2O temperature_k)
string(JSON thermopack_phase_kernel_state_count LENGTH "${thermopack_phase_kernel_json}" states)
if(NOT thermopack_phase_kernel_schema STREQUAL
       "MPMC_HNU/CPA/ThermoPack-phase-kernel-reference/v2")
    message(FATAL_ERROR "unexpected frozen ThermoPack CPA phase-kernel schema")
endif()
if(NOT thermopack_phase_kernel_commit STREQUAL
       "d68c794c7342bfc6938eb424a1fbb88b7780b738")
    message(FATAL_ERROR "unexpected ThermoPack phase-kernel revision")
endif()
if(NOT thermopack_phase_kernel_state_count EQUAL 5)
    message(FATAL_ERROR "ThermoPack CPA phase-kernel reference must contain five states")
endif()

function(read_phase_reference json index phase output_var)
    string(JSON phase_meoh GET "${json}" states ${index} ${phase} composition MEOH)
    string(JSON phase_rho GET "${json}" states ${index} ${phase} molar_density_mol_per_m3)
    string(JSON phase_z GET "${json}" states ${index} ${phase} compressibility_factor)
    string(JSON phase_lnphi_meoh GET "${json}" states ${index} ${phase} ln_phi MEOH)
    string(JSON phase_lnphi_h2o GET "${json}" states ${index} ${phase} ln_phi H2O)
    string(JSON phase_p_physical GET "${json}" states ${index} ${phase} pressure_physical_pa)
    string(JSON phase_p_association GET "${json}" states ${index} ${phase} pressure_association_pa)
    string(JSON phase_mu_cubic_meoh GET "${json}" states ${index} ${phase} mu_cubic_over_rt MEOH)
    string(JSON phase_mu_cubic_h2o GET "${json}" states ${index} ${phase} mu_cubic_over_rt H2O)
    string(JSON phase_mu_assoc_meoh GET "${json}" states ${index} ${phase} mu_association_over_rt MEOH)
    string(JSON phase_mu_assoc_h2o GET "${json}" states ${index} ${phase} mu_association_over_rt H2O)
    set(${output_var}
        "PhaseReference{${phase_meoh}, ${phase_rho}, ${phase_z}, ${phase_lnphi_meoh}, ${phase_lnphi_h2o}, ${phase_p_physical}, ${phase_p_association}, ${phase_mu_cubic_meoh}, ${phase_mu_cubic_h2o}, ${phase_mu_assoc_meoh}, ${phase_mu_assoc_h2o}}"
        PARENT_SCOPE)
endfunction()

set(thermopack_phase_kernel_initializers "")
foreach(index RANGE 0 4)
    string(JSON state_temperature GET "${thermopack_phase_kernel_json}" states ${index} temperature_k)
    string(JSON state_pressure GET "${thermopack_phase_kernel_json}" states ${index} pressure_pa)
    read_phase_reference("${thermopack_phase_kernel_json}" ${index} liquid liquid_initializer)
    read_phase_reference("${thermopack_phase_kernel_json}" ${index} vapor vapor_initializer)
    string(APPEND thermopack_phase_kernel_initializers
        "    ThermoPackPhaseKernelState{${state_temperature}, ${state_pressure}, ${liquid_initializer}, ${vapor_initializer}},\n")
endforeach()

set(thermopack_phase_kernel_generated_header
    "${CMAKE_CURRENT_BINARY_DIR}/cpa_thermopack_phase_kernel_generated.hpp")
file(GENERATE OUTPUT "${thermopack_phase_kernel_generated_header}" CONTENT
"#ifndef MPMC_TEST_CPA_THERMOPACK_PHASE_KERNEL_GENERATED_HPP\n#define MPMC_TEST_CPA_THERMOPACK_PHASE_KERNEL_GENERATED_HPP\n\n#include <array>\n\nnamespace cpa_thermopack_phase_kernel {\nstruct PhaseReference {\n    double composition_methanol;\n    double molar_density_mol_per_m3;\n    double compressibility_factor;\n    double ln_phi_methanol;\n    double ln_phi_water;\n    double pressure_physical_pa;\n    double pressure_association_pa;\n    double mu_cubic_over_rt_methanol;\n    double mu_cubic_over_rt_water;\n    double mu_association_over_rt_methanol;\n    double mu_association_over_rt_water;\n};\n\nstruct ThermoPackPhaseKernelState {\n    double temperature_k;\n    double pressure_pa;\n    PhaseReference liquid;\n    PhaseReference vapor;\n};\n\ninline constexpr double gas_constant_j_per_mol_k = ${thermopack_phase_kernel_rgas};\ninline constexpr double methanol_alpha_critical_temperature_k = ${thermopack_phase_kernel_tc_meoh};\ninline constexpr double water_alpha_critical_temperature_k = ${thermopack_phase_kernel_tc_h2o};\ninline constexpr std::array<ThermoPackPhaseKernelState, ${thermopack_phase_kernel_state_count}> states{{\n${thermopack_phase_kernel_initializers}}};\n\n} // namespace cpa_thermopack_phase_kernel\n\n#endif // MPMC_TEST_CPA_THERMOPACK_PHASE_KERNEL_GENERATED_HPP\n")


set(clapeyron_gate_e_file
    "${CMAKE_CURRENT_LIST_DIR}/../cpa_clapeyron_oracle/clapeyron_229b094_meoh_h2o_phase_kernel_33315k.json")
file(READ "${clapeyron_gate_e_file}" clapeyron_gate_e_json)
string(JSON clapeyron_gate_e_schema GET "${clapeyron_gate_e_json}" schema)
string(JSON clapeyron_gate_e_commit GET "${clapeyron_gate_e_json}" source commit)
string(JSON clapeyron_gate_e_profile GET "${clapeyron_gate_e_json}" model profile)
string(JSON clapeyron_gate_e_status GET "${clapeyron_gate_e_json}" model audit status)
string(JSON clapeyron_gate_e_rgas GET "${clapeyron_gate_e_json}" gas_constant_j_per_mol_k)
string(JSON clapeyron_gate_e_state_count LENGTH "${clapeyron_gate_e_json}" states)
if(NOT clapeyron_gate_e_schema STREQUAL
       "MPMC_HNU/CPA/Clapeyron-phase-kernel-reference/v1")
    message(FATAL_ERROR "unexpected frozen Clapeyron CPA oracle schema")
endif()
if(NOT clapeyron_gate_e_commit STREQUAL
       "229b09452f36c2f812486150df0bb43b197bb4e5")
    message(FATAL_ERROR "unexpected frozen Clapeyron revision")
endif()
if(NOT clapeyron_gate_e_profile STREQUAL
       "CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1" OR
   NOT clapeyron_gate_e_status STREQUAL "exactly_representable")
    message(FATAL_ERROR "frozen Clapeyron oracle no longer matches the audited CPA profile")
endif()
if(NOT clapeyron_gate_e_state_count EQUAL 10)
    message(FATAL_ERROR "Clapeyron Gate-E oracle must contain exactly ten states")
endif()

function(read_clapeyron_gate_e_state json index output_var)
    string(JSON state_temperature GET "${json}" states ${index} temperature_k)
    string(JSON state_rho GET "${json}" states ${index} molar_density_mol_per_m3)
    string(JSON state_meoh GET "${json}" states ${index} composition MEOH)
    string(JSON state_f_res GET "${json}" states ${index} f_res_reduced_extensive)
    string(JSON state_pressure GET "${json}" states ${index} pressure_pa)
    string(JSON state_mu_meoh GET "${json}" states ${index} mu_residual_over_rt MEOH)
    string(JSON state_mu_h2o GET "${json}" states ${index} mu_residual_over_rt H2O)
    string(JSON state_lnphi_meoh GET "${json}" states ${index} ln_phi MEOH)
    string(JSON state_lnphi_h2o GET "${json}" states ${index} ln_phi H2O)
    string(JSON state_assoc_residual
           GET "${json}" states ${index} association_equation_max_residual)
    set(${output_var}
        "ClapeyronState{${state_temperature}, ${state_rho}, ${state_meoh}, ${state_f_res}, ${state_pressure}, ${state_mu_meoh}, ${state_mu_h2o}, ${state_lnphi_meoh}, ${state_lnphi_h2o}, ${state_assoc_residual}}"
        PARENT_SCOPE)
endfunction()

set(clapeyron_gate_e_initializers "")
foreach(index RANGE 0 9)
    read_clapeyron_gate_e_state(
        "${clapeyron_gate_e_json}" ${index} clapeyron_state_initializer)
    string(APPEND clapeyron_gate_e_initializers
        "    ${clapeyron_state_initializer},\n")
endforeach()

set(clapeyron_gate_e_generated_header
    "${CMAKE_CURRENT_BINARY_DIR}/cpa_clapeyron_gate_e_generated.hpp")
file(GENERATE OUTPUT "${clapeyron_gate_e_generated_header}" CONTENT
"#ifndef MPMC_TEST_CPA_CLAPEYRON_GATE_E_GENERATED_HPP\n#define MPMC_TEST_CPA_CLAPEYRON_GATE_E_GENERATED_HPP\n\n#include <array>\n#include <string_view>\n\nnamespace cpa_clapeyron_gate_e {\nstruct ClapeyronState {\n    double temperature_k;\n    double molar_density_mol_per_m3;\n    double composition_methanol;\n    double f_res;\n    double pressure_pa;\n    double mu_res_over_rt_methanol;\n    double mu_res_over_rt_water;\n    double ln_phi_methanol;\n    double ln_phi_water;\n    double association_equation_max_residual;\n};\n\ninline constexpr std::string_view source_commit = \"${clapeyron_gate_e_commit}\";\ninline constexpr double gas_constant_j_per_mol_k = ${clapeyron_gate_e_rgas};\ninline constexpr std::array<ClapeyronState, ${clapeyron_gate_e_state_count}> states{{\n${clapeyron_gate_e_initializers}}};\n\n} // namespace cpa_clapeyron_gate_e\n\n#endif // MPMC_TEST_CPA_CLAPEYRON_GATE_E_GENERATED_HPP\n")

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

add_executable(mpmc_cpa_thermopack_phase_kernel_audit_tests
    cpa_thermopack_phase_kernel_audit_test.cpp)
target_link_libraries(mpmc_cpa_thermopack_phase_kernel_audit_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_cpa_thermopack_phase_kernel_audit_tests PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}")
set_target_properties(mpmc_cpa_thermopack_phase_kernel_audit_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_thermopack_parameter_identity_tests
    cpa_thermopack_parameter_identity_test.cpp)
target_link_libraries(mpmc_cpa_thermopack_parameter_identity_tests PRIVATE mpmc::flash)
target_include_directories(mpmc_cpa_thermopack_parameter_identity_tests PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}")
set_target_properties(mpmc_cpa_thermopack_parameter_identity_tests PROPERTIES CXX_EXTENSIONS OFF)


add_executable(mpmc_cpa_clapeyron_gate_e_tests
    cpa_clapeyron_gate_e_test.cpp)
target_link_libraries(mpmc_cpa_clapeyron_gate_e_tests PRIVATE
    mpmc::flash mpmc::ad)
target_include_directories(mpmc_cpa_clapeyron_gate_e_tests PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}")
set_target_properties(mpmc_cpa_clapeyron_gate_e_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_residual_helmholtz_ad_tests
    cpa_residual_helmholtz_ad_test.cpp)
target_link_libraries(mpmc_cpa_residual_helmholtz_ad_tests PRIVATE
    mpmc::flash mpmc::ad)
target_include_directories(mpmc_cpa_residual_helmholtz_ad_tests PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}")
set_target_properties(mpmc_cpa_residual_helmholtz_ad_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_residual_helmholtz_limits_tests
    cpa_residual_helmholtz_limits_test.cpp)
target_link_libraries(mpmc_cpa_residual_helmholtz_limits_tests PRIVATE
    mpmc::flash mpmc::ad)
target_include_directories(mpmc_cpa_residual_helmholtz_limits_tests PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}")
set_target_properties(mpmc_cpa_residual_helmholtz_limits_tests PROPERTIES CXX_EXTENSIONS OFF)

add_executable(mpmc_cpa_residual_helmholtz_temperature_tests
    cpa_residual_helmholtz_temperature_test.cpp)
target_link_libraries(mpmc_cpa_residual_helmholtz_temperature_tests PRIVATE
    mpmc::flash mpmc::ad)
target_include_directories(mpmc_cpa_residual_helmholtz_temperature_tests PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}")
set_target_properties(mpmc_cpa_residual_helmholtz_temperature_tests PROPERTIES CXX_EXTENSIONS OFF)

set(all_targets
    mpmc_cpa_physical_validation_tests
    mpmc_cpa_physical_validation_headers
    mpmc_cpa_thermopack_delta_tests
    mpmc_cpa_thermopack_phase_kernel_audit_tests
    mpmc_cpa_thermopack_parameter_identity_tests
    mpmc_cpa_clapeyron_gate_e_tests
    mpmc_cpa_residual_helmholtz_ad_tests
    mpmc_cpa_residual_helmholtz_limits_tests
    mpmc_cpa_residual_helmholtz_temperature_tests)

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

add_test(NAME flash.cpa_physical_validation.thermopack_phase_kernel_audit
         COMMAND mpmc_cpa_thermopack_phase_kernel_audit_tests)
set_tests_properties(flash.cpa_physical_validation.thermopack_phase_kernel_audit PROPERTIES
    LABELS "flash;cpa;physical_validation;thermopack;phase_kernel;diagnostic;full_physical"
    TIMEOUT 180)

add_test(NAME flash.cpa_physical_validation.thermopack_parameter_identity
         COMMAND mpmc_cpa_thermopack_parameter_identity_tests)
set_tests_properties(flash.cpa_physical_validation.thermopack_parameter_identity PROPERTIES
    LABELS "flash;cpa;physical_validation;thermopack;parameters;provenance;full_physical;sanitizer_safety"
    TIMEOUT 60)


add_test(NAME flash.cpa_physical_validation.clapeyron_gate_e
         COMMAND mpmc_cpa_clapeyron_gate_e_tests)
set_tests_properties(flash.cpa_physical_validation.clapeyron_gate_e PROPERTIES
    LABELS "flash;cpa;physical_validation;helmholtz;ad;clapeyron;external_oracle;full_physical;sanitizer_safety"
    TIMEOUT 60)

add_test(NAME flash.cpa_physical_validation.residual_helmholtz_ad
         COMMAND mpmc_cpa_residual_helmholtz_ad_tests)
set_tests_properties(flash.cpa_physical_validation.residual_helmholtz_ad PROPERTIES
    LABELS "flash;cpa;physical_validation;helmholtz;ad;thermopack;full_physical;sanitizer_safety"
    TIMEOUT 180)

add_test(NAME flash.cpa_physical_validation.residual_helmholtz_limits
         COMMAND mpmc_cpa_residual_helmholtz_limits_tests)
set_tests_properties(flash.cpa_physical_validation.residual_helmholtz_limits PROPERTIES
    LABELS "flash;cpa;physical_validation;helmholtz;ad;limits;full_physical;sanitizer_safety"
    TIMEOUT 60)

add_test(NAME flash.cpa_physical_validation.residual_helmholtz_temperature
         COMMAND mpmc_cpa_residual_helmholtz_temperature_tests)
set_tests_properties(flash.cpa_physical_validation.residual_helmholtz_temperature PROPERTIES
    LABELS "flash;cpa;physical_validation;helmholtz;ad;temperature;full_physical;sanitizer_safety"
    TIMEOUT 60)

add_test(NAME flash.cpa_physical_validation.headers
         COMMAND mpmc_cpa_physical_validation_headers)
set_tests_properties(flash.cpa_physical_validation.headers PROPERTIES
    LABELS "flash;cpa;physical_validation;headers;full_physical;sanitizer_safety"
    TIMEOUT 60)