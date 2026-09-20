cmake_minimum_required(VERSION 3.21)
project(mpmc_flow_discretization_petsc_tests LANGUAGES CXX)

add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/flow_discretization_petsc"
    flow-discretization-petsc-library)
enable_testing()

add_executable(
    mpmc_flow_discretization_petsc_tests
    distributed_component_conservation_test.cpp
    fixed_three_phase_snes_assembly_test.cpp
    single_phase_snes_assembly_test.cpp
    two_phase_snes_assembly_test.cpp
    variable_cardinality_natural_variable_numbering_test.cpp
    distributed_component_conservation_header.cpp
    distributed_energy_conservation_header.cpp
    energy_global_assembly_mapping_header.cpp
    complete_natural_variable_assembly_snapshot_header.cpp
    complete_natural_variable_petsc_materialization_header.cpp
    natural_variable_newton_linear_system_header.cpp
    natural_variable_snes_solver_header.cpp
    fixed_three_phase_snes_assembly_header.cpp
    single_phase_snes_assembly_header.cpp
    two_phase_snes_assembly_header.cpp
    variable_cardinality_natural_variable_numbering_header.cpp
    global_component_assembly_mapping_header.cpp
    fugacity_equilibrium_global_assembly_mapping_header.cpp)
target_link_libraries(
    mpmc_flow_discretization_petsc_tests
    PRIVATE mpmc::flow_discretization_petsc)
set_target_properties(
    mpmc_flow_discretization_petsc_tests
    PROPERTIES CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(
        mpmc_flow_discretization_petsc_tests
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
            -Werror)
endif()

add_test(
    NAME flow_discretization.petsc.distributed_owner_targeted_conservation
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:mpmc_flow_discretization_petsc_tests>
        ${MPIEXEC_POSTFLAGS})
set_tests_properties(
    flow_discretization.petsc.distributed_owner_targeted_conservation
    PROPERTIES
        LABELS "flow_discretization;flow_discretization.petsc;mpi;component_conservation;energy_conservation;global_assembly;complete_assembly;petsc_materialization;newton_linear_system;snes;gmres;asm;physical_assembly;owner_targeted"
        TIMEOUT 60
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1;OMPI_MCA_btl_vader_single_copy_mechanism=none")
