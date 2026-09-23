import json
import os
import re
import subprocess

router_path = ".github/workflows/pr_incremental_ci.yml"
all_ad_suites = ["arithmetic", "jacobian", "math", "runtime"]
ad_platforms = [
    {"name": "GCC Debug + ASan/UBSan", "os": "ubuntu-24.04", "runner": "ubuntu-24.04",
     "build_type": "Debug", "compiler": "-DCMAKE_CXX_COMPILER=g++", "sanitizer": "ON"},
    {"name": "Clang Release", "os": "ubuntu-24.04", "runner": "ubuntu-24.04",
     "build_type": "Release", "compiler": "-DCMAKE_CXX_COMPILER=clang++", "sanitizer": "OFF"},
    {"name": "MSVC Release", "os": "windows-2022", "runner": "windows-2022",
     "build_type": "Release", "compiler": "", "sanitizer": "OFF"},
]

all_thermo_suites = ["contracts", "pr76", "pr76_mixture", "pr76_pt"]
thermo_platforms = [
    {"name": "GCC Debug + ASan/UBSan", "os": "ubuntu-24.04", "runner": "ubuntu-24.04",
     "config": "Debug", "compiler": "-DCMAKE_CXX_COMPILER=g++", "sanitizer": "ON"},
    {"name": "Clang Release", "os": "ubuntu-24.04", "runner": "ubuntu-24.04",
     "config": "Release", "compiler": "-DCMAKE_CXX_COMPILER=clang++", "sanitizer": "OFF"},
    {"name": "MSVC Release", "os": "windows-2022", "runner": "windows-2022",
     "config": "Release", "compiler": "", "sanitizer": "OFF"},
]
thermo_targets = {
    "contracts": "mpmc_thermodynamics_contract_tests",
    "pr76": "mpmc_thermodynamics_pr76_tests",
    "pr76_mixture": "mpmc_thermodynamics_pr76_mixture_tests",
    "pr76_pt": "mpmc_thermodynamics_pr76_pt_tests",
}

ownership = {
    "flow_core": [
        "modules/flow/include/**",
        "modules/flow/CMakeLists.txt",
        "modules/thermodynamics/CMakeLists.txt",
        "modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_fugacity.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_density.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/pr76_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/sw92_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_*.hpp",
        "modules/ad/include/mpmc/ad/**",
        "modules/ad/CMakeLists.txt",
        "modules/flash/include/mpmc/flash/pt_phase_set.hpp",
        "modules/flash/include/mpmc/flash/pt_phase_transition.hpp",
        "modules/flash/include/mpmc/flash/pt_flash_backend.hpp",
        "tests/flow/core/**",
        ".github/workflows/flow_core.yml",
    ],
    "flow_discretization": [
        "modules/flow_discretization/**",
        "modules/flow/discretization/include/**",
        "modules/flow/discretization/CMakeLists.txt",
        "modules/well/include/**",
        "modules/well/CMakeLists.txt",
        "modules/well/discretization/include/**",
        "modules/well/discretization/CMakeLists.txt",
        "tests/well/discretization/**",
        "modules/flow/include/mpmc/flow/phase_potential_upwind.hpp",
        "modules/flow/include/mpmc/flow/phase_transport.hpp",
        "modules/flow/include/mpmc/flow/saturation_constitutive.hpp",
        "modules/flow/include/mpmc/flow/natural_variable_cell_state.hpp",
        "modules/flow/include/mpmc/flow/component_accumulation.hpp",
        "modules/flow/include/mpmc/flow/component_accumulation_time.hpp",
        "modules/flow/include/mpmc/flow/energy_accumulation.hpp",
        "modules/flow/CMakeLists.txt",
        "modules/discretization/include/mpmc/discretization/tpfa_internal_face_transmissibility_snapshot_3d.hpp",
        "modules/discretization/include/mpmc/discretization/tpfa_static_face_transmissibility_3d.hpp",
        "modules/discretization/include/mpmc/discretization/transmissibility_admissibility_3d.hpp",
        "modules/discretization/CMakeLists.txt",
        "modules/mesh/include/mpmc/mesh/entity.hpp",
        "modules/mesh/include/mpmc/mesh/partition_snapshot.hpp",
        "modules/mesh/include/mpmc/mesh/topology.hpp",
        "tests/flow_discretization/core/**",
        ".github/workflows/flow_discretization.yml",
    ],
    "flow_discretization_petsc": [
        "modules/flow_discretization_petsc/**",
        "modules/flow/discretization/petsc/**",
        "modules/well/discretization/petsc/**",
        "modules/flow/include/mpmc/flow/fugacity_equilibrium_residual.hpp",
        "modules/flow/include/mpmc/flow/fugacity_equilibrium_linearization.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_density.hpp",
        "modules/flow_discretization/include/mpmc/flow_discretization/local_energy_conservation_residual.hpp",
        "modules/flow/discretization/include/mpmc/flow_discretization/local_energy_conservation_residual.hpp",
        "modules/flow_discretization/include/mpmc/flow_discretization/energy_face_flux.hpp",
        "modules/flow/discretization/include/mpmc/flow_discretization/energy_face_flux.hpp",
        "modules/flow/include/mpmc/flow/energy_accumulation.hpp",
        "modules/flow_discretization/include/mpmc/flow_discretization/owned_multi_cell_component_conservation.hpp",
        "modules/flow/discretization/include/mpmc/flow_discretization/owned_multi_cell_component_conservation.hpp",
        "modules/flow_discretization/include/mpmc/flow_discretization/local_component_conservation_residual.hpp",
        "modules/flow/discretization/include/mpmc/flow_discretization/local_component_conservation_residual.hpp",
        "modules/flow_discretization/include/mpmc/flow_discretization/normalized_component_face_contribution.hpp",
        "modules/flow/discretization/include/mpmc/flow_discretization/normalized_component_face_contribution.hpp",
        "modules/discretization_petsc/**",
        "modules/discretization/petsc/**",
        "modules/mesh/include/mpmc/mesh/partition_snapshot.hpp",
        "tests/flow_discretization/petsc/**",
        "tests/flash/sw92_phase_assigned_pt/physical_sample6.hpp",
        "tests/support/sw92/test_support.hpp",
        ".github/workflows/flow_discretization_petsc.yml",
    ],
    "pt_flash_backend": [
        "modules/flash/include/mpmc/flash/pt_flash_backend.hpp",
        "modules/flash/include/mpmc/flash/pr76_pt_flash_backend.hpp",
        "modules/flash/include/mpmc/flash/sw92_profile_c_pt_flash_backend.hpp",
        "modules/flash/include/mpmc/flash/cpa_pt_flash_backend.hpp",
        "modules/flash/include/mpmc/flash/pr76_*.hpp",
        "modules/flash/include/mpmc/flash/pt_*.hpp",
        "modules/flash/include/mpmc/flash/sw92_*.hpp",
        "modules/flash/include/mpmc/flash/cpa_*.hpp",
        "modules/flash/include/mpmc/flash/detail/**",
        "modules/thermodynamics/include/mpmc/thermodynamics/**",
        "tests/flash/pt_flash_backend/**",
        "tests/flash/pt_split/boundary_references.hpp",
        "tests/flash/pt_split/dew_limit_references.hpp",
        "tests/flash/pt_split/reference_boundary_decimal.py",
        "tests/flash/pt_split/reference_dew_limit_decimal.py",
        "tests/flash/pr76_three_phase/synthetic_fixture.hpp",
        "tests/flash/cpa_max3/test_support.hpp",
        "tests/flash/cpa_physical_validation/**",
        "tests/flash/sw92_phase_assigned_pt/physical_sample6.hpp",
        "tests/support/sw92/test_support.hpp",
        ".github/workflows/pt_flash_backend.yml",
    ],
    "pt_stability": [
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/pr76_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_stability.hpp",
        "modules/flash/**/CMakeLists.txt",
        "tests/flash/stability/**",
        "modules/thermodynamics/**/*.hpp",
        "modules/thermodynamics/**/CMakeLists.txt",
        "modules/ad/include/mpmc/ad/**",
        "modules/ad/CMakeLists.txt",
        ".github/workflows/stability.yml",
    ],
    "pt_split": [
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/pr76_stability.hpp",
        "modules/flash/include/mpmc/flash/rachford_rice.hpp",
        "modules/flash/include/mpmc/flash/pt_split.hpp",
        "modules/flash/include/mpmc/flash/pr76_split.hpp",
        "modules/flash/include/mpmc/flash/pt_phase_set.hpp",
        "modules/flash/include/mpmc/flash/pt_vle_phase_set.hpp",
        "modules/flash/include/mpmc/flash/pr76_phase_set.hpp",
        "modules/flash/include/mpmc/flash/pt_sensitivity.hpp",
        "modules/flash/include/mpmc/flash/pr76_sensitivity.hpp",
        "modules/flash/include/mpmc/flash/detail/**",
        "modules/flash/CMakeLists.txt",
        "modules/flash/sensitivity/CMakeLists.txt",
        "modules/thermodynamics/**/*.hpp",
        "modules/thermodynamics/**/CMakeLists.txt",
        "modules/ad/**/*.hpp",
        "modules/ad/**/CMakeLists.txt",
        "tests/flash/pt_split/**",
        ".github/workflows/pt_split.yml",
    ],
    "pr76_pt_continuation": [
        "modules/flash/include/mpmc/flash/pr76_pt_continuation.hpp",
        "modules/flash/include/mpmc/flash/pr76_three_phase.hpp",
        "modules/flash/include/mpmc/flash/pr76_split.hpp",
        "modules/flash/include/mpmc/flash/pr76_stability.hpp",
        "modules/flash/include/mpmc/flash/pt_three_phase.hpp",
        "modules/flash/include/mpmc/flash/pt_split.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/pr*.hpp",
        "tests/flash/pr76_pt_continuation/**",
        "tests/flash/pr76_three_phase/synthetic_fixture.hpp",
        ".github/workflows/pr76_pt_continuation.yml",
    ],
    "pr76_max3": [
        "modules/flash/include/mpmc/flash/pt_three_phase.hpp",
        "modules/flash/include/mpmc/flash/pr76_three_phase.hpp",
        "modules/flash/include/mpmc/flash/pr76_max3_phase_set.hpp",
        "modules/flash/include/mpmc/flash/pr76_pt_flash_backend.hpp",
        "modules/flash/include/mpmc/flash/pr76_*.hpp",
        "modules/flash/include/mpmc/flash/pt_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/pr*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_fugacity.hpp",
        "modules/flow/CMakeLists.txt",
        "modules/flow/include/mpmc/flow/natural_variable_cell_state.hpp",
        "modules/flow/include/mpmc/flow/fugacity_equilibrium_residual.hpp",
        "modules/flow/include/mpmc/flow/thermodynamics_fugacity_adapters.hpp",
        "tests/support/flow/fugacity_adapter_regression.hpp",
        "tests/flash/pr76_three_phase/**",
        ".github/workflows/pr76_three_phase.yml",
    ],
    "cpa_stability": [
        "modules/flash/include/mpmc/flash/cpa_stability.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_*.hpp",
        "tests/flash/cpa_stability/**",
        ".github/workflows/cpa_stability.yml",
        ".github/workflows/_cpp_ctest.yml",
        "tests/support/cmake/cpp_test_options.cmake",
    ],
    "cpa_split": [
        "modules/flash/include/mpmc/flash/cpa_*.hpp",
        "modules/flash/include/mpmc/flash/pt_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_*.hpp",
        "tests/flash/cpa_split/**",
        "tests/flash/cpa_stability/**",
        "tests/thermodynamics/cpa_pt_phase/**",
        ".github/workflows/cpa_split.yml",
    ],
    "cpa_max3": [
        "modules/flash/include/mpmc/flash/cpa_*.hpp",
        "modules/flash/include/mpmc/flash/pt_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_fugacity.hpp",
        "modules/flow/CMakeLists.txt",
        "modules/flow/include/mpmc/flow/natural_variable_cell_state.hpp",
        "modules/flow/include/mpmc/flow/fugacity_equilibrium_residual.hpp",
        "modules/flow/include/mpmc/flow/thermodynamics_fugacity_adapters.hpp",
        "tests/support/flow/fugacity_adapter_regression.hpp",
        "tests/flash/cpa_max3/**",
        "tests/flash/cpa_split/**",
        ".github/workflows/cpa_max3.yml",
    ],
    "cpa_baseline": [
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_*.hpp",
        "tests/thermodynamics/cpa_baseline/**",
        ".github/workflows/cpa_baseline.yml",
        ".github/workflows/_cpp_ctest.yml",
        "tests/support/cmake/cpp_test_options.cmake",
    ],
    "cpa_pt_phase": [
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_fugacity.hpp",
        "modules/thermodynamics/cpa*.md",
        "tests/thermodynamics/cpa_pt_phase/**",
        ".github/workflows/cpa_pt_phase.yml",
        ".github/workflows/_cpp_ctest.yml",
        "tests/support/cmake/cpp_test_options.cmake",
    ],
    "cpa_physical_validation": [
        "tests/flash/cpa_thermopack_oracle/thermopack_d68c794_meoh_h2o_33315k.json",
        "tests/flash/cpa_thermopack_oracle/thermopack_d68c794_meoh_h2o_phase_kernel_33315k.json",
        "tests/flash/cpa_clapeyron_oracle/clapeyron_229b094_meoh_h2o_phase_kernel_33315k.json",
        "modules/thermodynamics/include/mpmc/thermodynamics/cpa_*.hpp",
        "modules/flash/include/mpmc/flash/cpa_*.hpp",
        "modules/flash/include/mpmc/flash/pt_*.hpp",
        "tests/flash/cpa_physical_validation/**",
        ".github/workflows/cpa_physical_validation.yml",
    ],
    "thermodynamics_contracts": [
        "modules/thermodynamics/**/*.hpp",
        "modules/thermodynamics/**/CMakeLists.txt",
        "tests/thermodynamics/**/*.cpp",
        "tests/thermodynamics/**/*.hpp",
        "tests/thermodynamics/**/*.py",
        "tests/thermodynamics/**/CMakeLists.txt",
        "modules/ad/include/mpmc/ad/**",
        "modules/ad/CMakeLists.txt",
        ".github/workflows/thermodynamics.yml",
    ],
    "physics_closure": [
        "modules/physics/**",
        "modules/flash/include/mpmc/flash/pr76_split.hpp",
        "modules/flash/include/mpmc/flash/pr76_sensitivity.hpp",
        "modules/flash/include/mpmc/flash/pt_sensitivity.hpp",
        "modules/flash/include/mpmc/flash/detail/**",
        "modules/flash/sensitivity/CMakeLists.txt",
        "modules/thermodynamics/**/*.hpp",
        "modules/thermodynamics/**/CMakeLists.txt",
        "modules/ad/**/*.hpp",
        "modules/ad/**/CMakeLists.txt",
        "tests/flash/pt_split/reference_sensitivity_decimal.py",
        "tests/physics/thermodynamic_closure/**",
        ".github/workflows/physics_closure.yml",
    ],
    "model_configuration": [
        "modules/model_configuration/**",
        "modules/thermodynamics/include/mpmc/thermodynamics/components.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/pr_parameters.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/pr76_*.hpp",
        "modules/thermodynamics/CMakeLists.txt",
        "modules/flash/CMakeLists.txt",
        "modules/flash/include/mpmc/flash/pr76_*.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/pt_split.hpp",
        "modules/flash/include/mpmc/flash/pt_three_phase.hpp",
        "modules/flash/include/mpmc/flash/rachford_rice.hpp",
        "modules/flash/include/mpmc/flash/pt_flash_backend.hpp",
        "modules/flash/include/mpmc/flash/pt_phase_set.hpp",
        "modules/flash/include/mpmc/flash/pt_phase_transition.hpp",
        "modules/flash/include/mpmc/flash/pt_vle_phase_set.hpp",
        "tests/model_configuration/parameters/**",
        "tests/model_configuration/solver_settings/**",
        "tests/model_configuration/executable_model/**",
        "tests/model_configuration/registry/**",
        "tests/flash/pr76_three_phase/synthetic_fixture.hpp",
        "tests/flash/pr76_three_phase/sour_gas_fixture.hpp",
        "tests/flash/pr76_three_phase/sour_gas_references.hpp",
        "tests/thermodynamics/contracts/**",
        ".github/workflows/model_configuration.yml",
    ],
    "sw92_thermodynamics": [
        "tests/support/sw92/test_support.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/sw92_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_fugacity.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/selected_phase_density.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/components.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/pr76_roots.hpp",
        "modules/thermodynamics/CMakeLists.txt",
        "modules/thermodynamics/sw92.md",
        "tests/thermodynamics/sw92/**",
        ".github/workflows/sw92.yml",
    ],
    "sw92_family_vle": [
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_stability.hpp",
        "modules/flash/include/mpmc/flash/rachford_rice.hpp",
        "modules/flash/include/mpmc/flash/pt_split.hpp",
        "modules/flash/include/mpmc/flash/sw92_split.hpp",
        "modules/flash/CMakeLists.txt",
        "modules/thermodynamics/**/*.hpp",
        "modules/thermodynamics/**/CMakeLists.txt",
        "tests/flash/sw92_family_vle/**",
        ".github/workflows/sw92_family_vle.yml",
    ],
    "sw92_profile_c_pt": [
        "tests/support/sw92/test_support.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_split.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_no_w.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_pt.hpp",
        "modules/flash/sw92_phase_assigned_pt.md",
        "modules/flash/sw92_authoritative_three_phase_audit.md",
        "modules/thermodynamics/**/*.hpp",
        "tests/flash/sw92_family_vle/**",
        "tests/flash/sw92_phase_assigned_joint/**",
        "tests/flash/sw92_phase_assigned_h_side_witness/**",
        "tests/flash/sw92_phase_assigned_three_phase/**",
        "tests/flash/sw92_phase_assigned_three_phase_closure/**",
        "tests/flash/sw92_phase_assigned_no_w/**",
        "tests/flash/sw92_phase_assigned_pt/**",
        ".github/workflows/sw92_phase_assigned_pt.yml",
    ],              "sw92_profile_c_phase_set": [
        "tests/support/sw92/test_support.hpp",
        "modules/flash/include/mpmc/flash/sw92_profile_c_phase_set.hpp",
        "modules/flash/include/mpmc/flash/sw92_pt_flash.hpp",
        "modules/flash/include/mpmc/flash/sw92_profile_c_pt_flash_backend.hpp",
        "modules/flash/sw92_profile_c_phase_set.md",
        "modules/flash/sw92_pt_flash_correctness.md",
        "modules/flash/include/mpmc/flash/pt_phase_set.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_boundary.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_pt.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_no_w.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp",
        "modules/thermodynamics/**/*.hpp",
        "modules/flow/CMakeLists.txt",
        "modules/flow/include/mpmc/flow/natural_variable_cell_state.hpp",
        "modules/flow/include/mpmc/flow/fugacity_equilibrium_residual.hpp",
        "modules/flow/include/mpmc/flow/thermodynamics_fugacity_adapters.hpp",
        "tests/support/flow/fugacity_adapter_regression.hpp",
        "tests/flash/sw92_profile_c_phase_set/**",
        "tests/flash/sw92_phase_assigned_boundary/**",
        "tests/flash/sw92_phase_assigned_pt/**",
        "tests/flash/pt_split/phase_set_test.cpp",
        "tests/flash/pt_split/CMakeLists.txt",
        "README.md",
        ".github/workflows/sw92_profile_c_phase_set.yml",
    ],
    "sw92_profile_c_sensitivity": [
        "tests/support/sw92/test_support.hpp",
        "modules/ad/**",
        "modules/thermodynamics/include/mpmc/thermodynamics/sw92_*.hpp",
        "modules/flash/include/mpmc/flash/pt_sensitivity.hpp",
        "modules/flash/include/mpmc/flash/detail/pt_sensitivity_detail.hpp",
        "modules/flash/include/mpmc/flash/detail/sw92_profile_c_sensitivity_local.hpp",
        "modules/flash/include/mpmc/flash/sw92_profile_c_sensitivity.hpp",
        "modules/flash/include/mpmc/flash/sw92_profile_c_phase_set.hpp",
        "tests/flash/sw92_profile_c_sensitivity/**",
        "tests/flash/sw92_phase_assigned_pt/physical_sample6.hpp",
        "tests/flash/sw92_phase_assigned_pt/reference_physical_sample6_decimal.py",
        "tests/thermodynamics/sw92/**",
        ".github/workflows/sw92_profile_c_sensitivity.yml",
    ],
    "sw92_physics_closure": [
        "modules/ad/**",
        "modules/physics/include/mpmc/physics/thermodynamic_closure.hpp",
        "modules/physics/include/mpmc/physics/sw92_thermodynamic_closure.hpp",
        "modules/physics/CMakeLists.txt",
        "modules/flash/include/mpmc/flash/pt_sensitivity.hpp",
        "modules/flash/include/mpmc/flash/detail/pt_sensitivity_detail.hpp",
        "modules/flash/include/mpmc/flash/detail/sw92_profile_c_sensitivity_local.hpp",
        "modules/flash/include/mpmc/flash/sw92_profile_c_sensitivity.hpp",
        "modules/flash/include/mpmc/flash/sw92_profile_c_phase_set.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_*.hpp",
        "modules/thermodynamics/include/mpmc/thermodynamics/sw92_*.hpp",
        "tests/physics/sw92_thermodynamic_closure/**",
        "tests/flash/sw92_profile_c_sensitivity/**",
        "tests/flash/sw92_phase_assigned_pt/physical_sample6.hpp",
        "tests/flash/sw92_phase_assigned_pt/reference_physical_sample6_decimal.py",
        "tests/support/sw92/test_support.hpp",
        ".github/workflows/physics_sw92_closure.yml",
    ],
    "sw92_phase_assigned_joint": [
        "tests/support/sw92/test_support.hpp",
        "modules/flash/include/mpmc/flash/rachford_rice.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp",
        "modules/flash/sw92_phase_assigned_joint.md",
        "modules/flash/sw92_phase_assigned_c2_stability_audit.md",
        "modules/thermodynamics/**/*.hpp",
        "tests/flash/sw92_phase_assigned_joint/**",
        ".github/workflows/sw92_phase_assigned_joint.yml",
    ],
    "sw92_phase_assigned_h_side_witness": [
        "tests/support/sw92/test_support.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp",
        "modules/flash/sw92_phase_assigned_h_side_witness.md",
        "modules/thermodynamics/**/*.hpp",
        "tests/flash/sw92_phase_assigned_joint/**",
        "tests/flash/sw92_phase_assigned_h_side_witness/**",
        ".github/workflows/sw92_phase_assigned_h_side_witness.yml",
    ],
    "sw92_phase_assigned_no_w": [
        "tests/support/sw92/test_support.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/pt_split.hpp",
        "modules/flash/include/mpmc/flash/sw92_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_split.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_no_w.hpp",
        "modules/flash/sw92_phase_assigned_no_w.md",
        "modules/thermodynamics/**/*.hpp",
        "tests/flash/sw92_family_vle/**",
        "tests/flash/sw92_phase_assigned_joint/**",
        "tests/flash/sw92_phase_assigned_no_w/**",
        ".github/workflows/sw92_phase_assigned_no_w.yml",
    ],
    "sw92_phase_assigned_three_phase": [
        "tests/support/sw92/test_support.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase.hpp",
        "modules/flash/sw92_phase_assigned_three_phase.md",
        "modules/thermodynamics/**/*.hpp",
        "tests/flash/sw92_asymmetric_max2_physical_ternary/reference_decimal.py",
        "tests/flash/sw92_phase_assigned_joint/**",
        "tests/flash/sw92_phase_assigned_h_side_witness/**",
        "tests/flash/sw92_phase_assigned_three_phase/**",
        ".github/workflows/sw92_phase_assigned_three_phase.yml",
        ".github/workflows/_sw92_phase_assigned_topology.yml",
    ],
    "sw92_phase_assigned_three_phase_closure": [
        "tests/support/sw92/test_support.hpp",
        "modules/flash/include/mpmc/flash/pt_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_stability.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp",
        "modules/flash/sw92_phase_assigned_three_phase_closure.md",
        "modules/thermodynamics/**/*.hpp",
        "tests/flash/sw92_phase_assigned_joint/**",
        "tests/flash/sw92_phase_assigned_h_side_witness/**",
        "tests/flash/sw92_phase_assigned_three_phase/**",
        "tests/flash/sw92_phase_assigned_three_phase_closure/**",
        ".github/workflows/sw92_phase_assigned_three_phase_closure.yml",
        ".github/workflows/_sw92_phase_assigned_topology.yml",
    ],
    "sw92_phase_assigned_boundary": [
        "tests/support/sw92/test_support.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_boundary.hpp",
        "modules/flash/sw92_phase_assigned_boundary.md",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_pt.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_no_w.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp",
        "modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp",
        "modules/thermodynamics/**/*.hpp",
        "tests/flash/sw92_phase_assigned_boundary/**",
        "tests/flash/sw92_phase_assigned_pt/**",
        "tests/flash/sw92_phase_assigned_three_phase_closure/**",
        "tests/flash/sw92_phase_assigned_no_w/**",
        "tests/flash/sw92_phase_assigned_joint/**",
        "tests/flash/sw92_phase_assigned_h_side_witness/**",
        ".github/workflows/sw92_phase_assigned_boundary.yml",
        ".github/workflows/_sw92_phase_assigned_topology.yml",
    ],

}

def compile_pattern(pattern):
    out = ["^"]
    i = 0
    while i < len(pattern):
        if pattern[i:i + 3] == "**/":
            out.append("(?:.*/)?")
            i += 3
        elif pattern[i:i + 2] == "**":
            out.append(".*")
            i += 2
        elif pattern[i] == "*":
            out.append("[^/]*")
            i += 1
        elif pattern[i] == "?":
            out.append("[^/]")
            i += 1
        else:
            out.append(re.escape(pattern[i]))
            i += 1
    out.append("$")
    return re.compile("".join(out))

compiled = {
    gate: [compile_pattern(pattern) for pattern in patterns]
    for gate, patterns in ownership.items()
}

def route(paths):
    result = {
        gate: any(
            pattern.match(path)
            for path in paths
            for pattern in patterns
        )
        for gate, patterns in compiled.items()
    }
    return result

def expected_route(**overrides):
    expected = {gate: False for gate in ownership}
    expected.update(overrides)
    return expected

def ad_suites_for(paths):
    selected = set()
    infrastructure = False
    shared = {
        "modules/ad/include/mpmc/ad/dual.hpp",
        "modules/ad/CMakeLists.txt",
        "CMakeLists.txt",
        "CMakePresets.json",
    }
    owned_exact = {
        "CMakeLists.txt",
        "CMakePresets.json",
        ".github/workflows/ad.yml",
    }
    for path in paths:
        owned = (
            path.startswith("modules/ad/") or
            path.startswith("tests/ad/") or
            path in owned_exact
        )
        if not owned:
            continue
        if path.endswith(".md"):
            continue
        if path == ".github/workflows/ad.yml":
            infrastructure = True
        elif path in shared:
            selected.update(all_ad_suites)
        elif (
            path == "modules/ad/include/mpmc/ad/runtime_differentiate.hpp" or
            path.startswith("tests/ad/runtime/")
        ):
            selected.add("runtime")
        elif path == "modules/ad/include/mpmc/ad/differentiate.hpp":
            selected.update(("jacobian", "runtime"))
        elif path.startswith("tests/ad/jacobian/"):
            selected.add("jacobian")
        elif path == "modules/ad/include/mpmc/ad/math.hpp":
            selected.update(("jacobian", "math", "runtime"))
        elif path.startswith("tests/ad/math/"):
            selected.add("math")
        elif path.startswith("tests/ad/"):
            selected.add("arithmetic")
        else:
            selected.update(all_ad_suites)
    if infrastructure and not selected:
        selected.update(all_ad_suites)
    return sorted(selected)

def ad_matrix(selected):
    return {
        "include": [
            {**platform, "suite": suite}
            for suite in selected
            for platform in ad_platforms
        ]
    }

thermo_prefix = "modules/thermodynamics/include/mpmc/thermodynamics/"

def thermo_suites_for(paths):
    selected = set()
    infrastructure = False
    for path in paths:
        if path.endswith(".md"):
            continue
        if path == ".github/workflows/thermodynamics.yml":
            infrastructure = True
        elif path == thermo_prefix + "selected_phase_fugacity.hpp":
            selected.add("pr76_pt")
        elif (
            path in [
                thermo_prefix + "pr76_phase.hpp",
                thermo_prefix + "pr76_roots.hpp",
            ]
            or path.startswith("tests/thermodynamics/pr76_pt/")
        ):
            selected.add("pr76_pt")
        elif path == thermo_prefix + "pr76_mixture.hpp":
            selected.update(("pr76_mixture", "pr76_pt"))
        elif path.startswith("tests/thermodynamics/pr76_mixture/"):
            selected.add("pr76_mixture")
        elif path == thermo_prefix + "pr76_pure.hpp":
            selected.update(("pr76", "pr76_mixture", "pr76_pt"))
        elif path.startswith("tests/thermodynamics/pr76/"):
            selected.add("pr76")
        elif path.startswith("tests/thermodynamics/contracts/"):
            selected.add("contracts")
        elif path.startswith("modules/ad/"):
            selected.update(("pr76", "pr76_mixture", "pr76_pt"))
        else:
            selected.update(all_thermo_suites)
    if infrastructure and not selected:
        selected.update(all_thermo_suites)
    return sorted(selected)

def thermo_matrix(selected):
    return {
        "include": [
            {
                **platform,
                "suite": suite,
                "target": thermo_targets[suite],
            }
            for suite in selected
            for platform in thermo_platforms
        ]
    }

checks = [
    (["tests/flow_discretization/petsc/scanner.cpp"], expected_route(flow_discretization_petsc=True)),
    (["modules/flow/discretization/include/mpmc/flow_discretization/cell_source.hpp"], expected_route(flow_discretization=True)),
    (["modules/well/discretization/include/mpmc/well_discretization/hydraulic_conductance.hpp"], expected_route(flow_discretization=True)),
    (["tests/well/discretization/hydraulic_conductance_test.cpp"], expected_route(flow_discretization=True)),
    (["modules/flow/discretization/petsc/include/mpmc/flow_discretization_petsc/physical_timestep_driver.hpp"], expected_route(flow_discretization_petsc=True)),
    (["modules/flow_discretization/include/mpmc/flow_discretization/cell_source.hpp"], expected_route(flow_discretization=True)),
    (["modules/flow_discretization_petsc/include/mpmc/flow_discretization_petsc/physical_timestep_driver.hpp"], expected_route(flow_discretization_petsc=True)),
    (["modules/discretization/petsc/include/mpmc/discretization_petsc/adapter.hpp"], expected_route(flow_discretization_petsc=True)),
    (["modules/flow/include/mpmc/flow/phase_transport.hpp"], expected_route(flow_core=True, flow_discretization=True)),
    (["modules/flash/include/mpmc/flash/pt_flash_backend.hpp"], expected_route(flow_core=True, pt_flash_backend=True, pr76_max3=True, cpa_split=True, cpa_max3=True, cpa_physical_validation=True, model_configuration=True)),
    (["tests/flash/stability/stability_test.cpp"], expected_route(pt_stability=True)),
    (["modules/flash/include/mpmc/flash/rachford_rice.hpp"], expected_route(pt_split=True, model_configuration=True, sw92_family_vle=True, sw92_phase_assigned_joint=True)),
    (["modules/thermodynamics/include/mpmc/thermodynamics/pr76.hpp"], expected_route(pt_flash_backend=True, pt_stability=True, pt_split=True, pr76_pt_continuation=True, pr76_max3=True, thermodynamics_contracts=True, physics_closure=True, sw92_family_vle=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_phase_assigned_joint=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_no_w=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["tests/flash/pr76_pt_continuation/continuation_test.cpp"], expected_route(pr76_pt_continuation=True)),
    (["tests/flash/pr76_three_phase/max3_test.cpp"], expected_route(pr76_max3=True)),
    (["modules/flash/include/mpmc/flash/pr76_pt_continuation.hpp"], expected_route(pt_flash_backend=True, pr76_pt_continuation=True, pr76_max3=True, model_configuration=True)),
    (["tests/flash/cpa_stability/stability_test.cpp"], expected_route(cpa_stability=True, cpa_split=True)),
    (["tests/flash/cpa_split/split_test.cpp"], expected_route(cpa_split=True, cpa_max3=True)),
    (["tests/flash/cpa_max3/max3_test.cpp"], expected_route(cpa_max3=True)),
    (["modules/flash/include/mpmc/flash/cpa_stability.hpp"], expected_route(pt_flash_backend=True, cpa_stability=True, cpa_split=True, cpa_max3=True, cpa_physical_validation=True)),
    ([".github/workflows/_cpp_ctest.yml"], expected_route(cpa_stability=True, cpa_baseline=True, cpa_pt_phase=True)),
    (["tests/thermodynamics/cpa_baseline/baseline_test.cpp"], expected_route(cpa_baseline=True, thermodynamics_contracts=True)),
    (["tests/thermodynamics/cpa_pt_phase/phase_test.cpp"], expected_route(cpa_split=True, cpa_pt_phase=True, thermodynamics_contracts=True)),
    (["tests/flash/cpa_physical_validation/physical_validation_test.cpp"], expected_route(pt_flash_backend=True, cpa_physical_validation=True)),
    (["modules/thermodynamics/include/mpmc/thermodynamics/cpa_phase.hpp"], expected_route(flow_core=True, pt_flash_backend=True, pt_stability=True, pt_split=True, cpa_stability=True, cpa_split=True, cpa_max3=True, cpa_baseline=True, cpa_pt_phase=True, cpa_physical_validation=True, thermodynamics_contracts=True, physics_closure=True, sw92_family_vle=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_phase_assigned_joint=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_no_w=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["modules/thermodynamics/include/mpmc/thermodynamics/components.hpp"], expected_route(pt_flash_backend=True, pt_stability=True, pt_split=True, thermodynamics_contracts=True, physics_closure=True, model_configuration=True, sw92_thermodynamics=True, sw92_family_vle=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_phase_assigned_joint=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_no_w=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["tests/thermodynamics/contracts/contracts_test.cpp"], expected_route(thermodynamics_contracts=True, model_configuration=True)),
    (["tests/thermodynamics/pr76_pt/pt_test.cpp"], expected_route(thermodynamics_contracts=True)),
    (["modules/physics/include/mpmc/physics/component_inventory.hpp"], expected_route(physics_closure=True)),
    (["tests/physics/thermodynamic_closure/closure_test.cpp"], expected_route(physics_closure=True)),
    (["modules/model_configuration/include/mpmc/model_configuration/model.hpp"], expected_route(model_configuration=True)),
    (["tests/model_configuration/parameters/parameter_test.cpp"], expected_route(model_configuration=True)),
    ([".github/workflows/thermodynamics.yml"], expected_route(thermodynamics_contracts=True)),
    ([".github/workflows/physics_closure.yml"], expected_route(physics_closure=True)),
    ([".github/workflows/model_configuration.yml"], expected_route(model_configuration=True)),
    (["modules/thermodynamics/include/mpmc/thermodynamics/sw92_phase.hpp"], expected_route(flow_core=True, pt_flash_backend=True, pt_stability=True, pt_split=True, thermodynamics_contracts=True, physics_closure=True, sw92_thermodynamics=True, sw92_family_vle=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_profile_c_sensitivity=True, sw92_physics_closure=True, sw92_phase_assigned_joint=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_no_w=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["tests/thermodynamics/sw92/sw92_test.cpp"], expected_route(thermodynamics_contracts=True, sw92_thermodynamics=True, sw92_profile_c_sensitivity=True)),
    (["tests/flash/sw92_family_vle/family_vle_test.cpp"], expected_route(sw92_family_vle=True, sw92_profile_c_pt=True, sw92_phase_assigned_no_w=True)),
    (["tests/flash/sw92_phase_assigned_pt/pt_test.cpp"], expected_route(sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_phase_assigned_boundary=True)),
    (["tests/support/sw92/test_support.hpp"], expected_route(flow_discretization_petsc=True, pt_flash_backend=True, sw92_thermodynamics=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_profile_c_sensitivity=True, sw92_physics_closure=True, sw92_phase_assigned_joint=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_no_w=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    ([".github/workflows/sw92.yml"], expected_route(sw92_thermodynamics=True)),
    ([".github/workflows/sw92_family_vle.yml"], expected_route(sw92_family_vle=True)),
    ([".github/workflows/sw92_phase_assigned_pt.yml"], expected_route(sw92_profile_c_pt=True)),
    ([".github/workflows/pr_incremental_ci.yml"], expected_route()),
    (["README.md"], expected_route(sw92_profile_c_phase_set=True)),
    (["modules/flash/include/mpmc/flash/sw92_profile_c_phase_set.hpp"], expected_route(pt_flash_backend=True, sw92_profile_c_phase_set=True, sw92_profile_c_sensitivity=True, sw92_physics_closure=True)),
    (["tests/flash/sw92_profile_c_phase_set/publication_test.cpp"], expected_route(sw92_profile_c_phase_set=True)),
    (["tests/flash/sw92_profile_c_sensitivity/sensitivity_test.cpp"], expected_route(sw92_profile_c_sensitivity=True, sw92_physics_closure=True)),
    (["modules/physics/include/mpmc/physics/sw92_thermodynamic_closure.hpp"], expected_route(physics_closure=True, sw92_physics_closure=True)),
    (["tests/physics/sw92_thermodynamic_closure/closure_test.cpp"], expected_route(sw92_physics_closure=True)),
    (["modules/ad/include/mpmc/ad/dual.hpp"], expected_route(flow_core=True, pt_stability=True, pt_split=True, thermodynamics_contracts=True, physics_closure=True, sw92_profile_c_sensitivity=True, sw92_physics_closure=True)),
    ([".github/workflows/sw92_profile_c_phase_set.yml"], expected_route(sw92_profile_c_phase_set=True)),
    ([".github/workflows/sw92_profile_c_sensitivity.yml"], expected_route(sw92_profile_c_sensitivity=True)),
    ([".github/workflows/physics_sw92_closure.yml"], expected_route(sw92_physics_closure=True)),
    (["modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp"], expected_route(pt_flash_backend=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_physics_closure=True, sw92_phase_assigned_joint=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["tests/flash/sw92_phase_assigned_joint/joint_test.cpp"], expected_route(sw92_profile_c_pt=True, sw92_phase_assigned_joint=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_no_w=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp"], expected_route(pt_flash_backend=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_physics_closure=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["tests/flash/sw92_phase_assigned_h_side_witness/witness_test.cpp"], expected_route(sw92_profile_c_pt=True, sw92_phase_assigned_h_side_witness=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["modules/flash/include/mpmc/flash/sw92_phase_assigned_no_w.hpp"], expected_route(pt_flash_backend=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_physics_closure=True, sw92_phase_assigned_no_w=True, sw92_phase_assigned_boundary=True)),
    (["tests/flash/sw92_phase_assigned_no_w/no_w_test.cpp"], expected_route(sw92_profile_c_pt=True, sw92_phase_assigned_no_w=True, sw92_phase_assigned_boundary=True)),
    ([".github/workflows/sw92_phase_assigned_joint.yml"], expected_route(sw92_phase_assigned_joint=True)),
    ([".github/workflows/sw92_phase_assigned_h_side_witness.yml"], expected_route(sw92_phase_assigned_h_side_witness=True)),
    ([".github/workflows/sw92_phase_assigned_no_w.yml"], expected_route(sw92_phase_assigned_no_w=True)),
    (["modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase.hpp"], expected_route(pt_flash_backend=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_physics_closure=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["tests/flash/sw92_phase_assigned_three_phase/three_phase_test.cpp"], expected_route(sw92_profile_c_pt=True, sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True)),
    (["modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp"], expected_route(pt_flash_backend=True, sw92_profile_c_pt=True, sw92_profile_c_phase_set=True, sw92_physics_closure=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["tests/flash/sw92_phase_assigned_three_phase_closure/closure_test.cpp"], expected_route(sw92_profile_c_pt=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
    (["modules/flash/include/mpmc/flash/sw92_phase_assigned_boundary.hpp"], expected_route(pt_flash_backend=True, sw92_profile_c_phase_set=True, sw92_physics_closure=True, sw92_phase_assigned_boundary=True)),
    (["tests/flash/sw92_phase_assigned_boundary/boundary_test.cpp"], expected_route(sw92_profile_c_phase_set=True, sw92_phase_assigned_boundary=True)),
    ([".github/workflows/sw92_phase_assigned_three_phase.yml"], expected_route(sw92_phase_assigned_three_phase=True)),
    ([".github/workflows/sw92_phase_assigned_three_phase_closure.yml"], expected_route(sw92_phase_assigned_three_phase_closure=True)),
    ([".github/workflows/sw92_phase_assigned_boundary.yml"], expected_route(sw92_phase_assigned_boundary=True)),
    ([".github/workflows/_sw92_phase_assigned_topology.yml"], expected_route(sw92_phase_assigned_three_phase=True, sw92_phase_assigned_three_phase_closure=True, sw92_phase_assigned_boundary=True)),
]
for paths, expected in checks:
    actual = route(paths)
    if actual != expected:
        raise RuntimeError(
            f"impact-router regression failed for {paths}: {actual} != {expected}"
        )
print(f"Impact-router regressions: {len(checks)}/{len(checks)} passed")

ad_checks = [
    (["modules/ad/include/mpmc/ad/math.hpp"], ["jacobian", "math", "runtime"]),
    (["tests/ad/math/math_test.cpp", ".github/workflows/ad.yml"], ["math"]),
    (["modules/ad/include/mpmc/ad/dual.hpp"], all_ad_suites),
    (["tests/ad/dual_test.cpp"], ["arithmetic"]),
    ([".github/workflows/ad.yml"], all_ad_suites),
    (["modules/ad/include/mpmc/ad/unknown.hpp"], all_ad_suites),
    (["README.md"], []),
    (["tests/flow_discretization/petsc/scanner.cpp"], []),
    (["modules/ad/include/mpmc/ad/differentiate.hpp"], ["jacobian", "runtime"]),
    (["tests/ad/jacobian/jacobian_test.cpp", ".github/workflows/ad.yml"], ["jacobian"]),
    (["modules/ad/include/mpmc/ad/differentiate.hpp", "AGENTS.md"], ["jacobian", "runtime"]),
    (["tests/ad/math/math_test.cpp", "tests/ad/jacobian/header_self_contained.cpp"],
     ["jacobian", "math"]),
    (["CMakeLists.txt"], all_ad_suites),
    (["modules/ad/include/mpmc/ad/runtime_differentiate.hpp"], ["runtime"]),
    (["tests/ad/runtime/runtime_test.cpp", ".github/workflows/ad.yml"], ["runtime"]),
    (["tests/ad/runtime/header_self_contained.cpp"], ["runtime"]),
    (["modules/ad/include/mpmc/ad/runtime_differentiate.hpp", "AGENTS.md"], ["runtime"]),
    (["modules/ad/include/mpmc/ad/math.hpp", "tests/ad/runtime/runtime_test.cpp"],
     ["jacobian", "math", "runtime"]),
    (["modules/ad/runtime_differentiate.md"], []),
]
for paths, expected in ad_checks:
    actual = ad_suites_for(paths)
    if actual != expected:
        raise RuntimeError(
            f"AD selector regression failed for {paths}: {actual} != {expected}"
        )
print(f"AD selector regressions: {len(ad_checks)}/{len(ad_checks)} passed")

thermo_checks = [
    ([thermo_prefix + "pr76_pure.hpp"], ["pr76", "pr76_mixture", "pr76_pt"]),
    (["tests/thermodynamics/pr76/pr76_test.cpp"], ["pr76"]),
    (["tests/thermodynamics/pr76/CMakeLists.txt"], ["pr76"]),
    ([thermo_prefix + "pr76_pure.hpp", ".github/workflows/thermodynamics.yml"],
     ["pr76", "pr76_mixture", "pr76_pt"]),
    ([thermo_prefix + "components.hpp"], all_thermo_suites),
    ([thermo_prefix + "pr_parameters.hpp"], all_thermo_suites),
    (["modules/thermodynamics/CMakeLists.txt"], all_thermo_suites),
    (["tests/thermodynamics/contracts/contracts_test.cpp"], ["contracts"]),
    (["modules/ad/include/mpmc/ad/dual.hpp"], ["pr76", "pr76_mixture", "pr76_pt"]),
    (["modules/ad/include/mpmc/ad/math.hpp"], ["pr76", "pr76_mixture", "pr76_pt"]),
    (["modules/ad/include/mpmc/ad/differentiate.hpp"], ["pr76", "pr76_mixture", "pr76_pt"]),
    (["modules/ad/include/mpmc/ad/runtime_differentiate.hpp"],
     ["pr76", "pr76_mixture", "pr76_pt"]),
    (["modules/ad/CMakeLists.txt"], ["pr76", "pr76_mixture", "pr76_pt"]),
    ([".github/workflows/thermodynamics.yml"], all_thermo_suites),
    (["modules/thermodynamics/README.md", "README.md"], []),
    ([thermo_prefix + "unknown.hpp"], all_thermo_suites),
    (["tests/thermodynamics/contracts/contracts_test.cpp", thermo_prefix + "pr76_pure.hpp"],
     all_thermo_suites),
    ([thermo_prefix + "pr76_mixture.hpp"], ["pr76_mixture", "pr76_pt"]),
    (["tests/thermodynamics/pr76_mixture/mixture_test.cpp"], ["pr76_mixture"]),
    (["tests/thermodynamics/pr76_mixture/CMakeLists.txt"], ["pr76_mixture"]),
    ([thermo_prefix + "pr76_mixture.hpp", ".github/workflows/thermodynamics.yml"],
     ["pr76_mixture", "pr76_pt"]),
    (["tests/thermodynamics/pr76/pr76_test.cpp", thermo_prefix + "pr76_mixture.hpp"],
     ["pr76", "pr76_mixture", "pr76_pt"]),
    (["modules/thermodynamics/pr76_mixture.md"], []),
    ([thermo_prefix + "pr76_phase.hpp"], ["pr76_pt"]),
    ([thermo_prefix + "selected_phase_fugacity.hpp"], ["pr76_pt"]),
    ([thermo_prefix + "pr76_roots.hpp"], ["pr76_pt"]),
    (["tests/thermodynamics/pr76_pt/pt_test.cpp"], ["pr76_pt"]),
    (["tests/thermodynamics/pr76_pt/test_support.hpp"], ["pr76_pt"]),
    (["tests/thermodynamics/pr76_pt/reference_decimal.py"], ["pr76_pt"]),
    (["tests/thermodynamics/pr76_pt/CMakeLists.txt"], ["pr76_pt"]),
    ([thermo_prefix + "pr76_phase.hpp", ".github/workflows/thermodynamics.yml"],
     ["pr76_pt"]),
    (["modules/thermodynamics/pr76_phase.md"], []),
]
for paths, expected in thermo_checks:
    actual = thermo_suites_for(paths)
    if actual != expected:
        raise RuntimeError(
            f"Thermodynamics selector regression failed for {paths}: {actual} != {expected}"
        )
print(
    f"Thermodynamics selector regressions: "
    f"{len(thermo_checks)}/{len(thermo_checks)} passed"
)

