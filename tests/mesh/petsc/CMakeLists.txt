cmake_minimum_required(VERSION 3.21)
project(mpmc_mesh_discretization_petsc_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/mesh" mesh-core)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/mesh_petsc" mesh-petsc)
add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/discretization_petsc"
    discretization-petsc)
add_subdirectory(
    "${CMAKE_CURRENT_LIST_DIR}/../../../modules/flow_discretization"
    flow-discretization-library)
enable_testing()

add_executable(mpmc_mesh_discretization_petsc_tests
    petsc_adapter_test.cpp
    adapter_header_self_contained.cpp
    discretization_adapter_header_self_contained.cpp
    flow_discretization_schedule_bridge_header.cpp)
target_link_libraries(
    mpmc_mesh_discretization_petsc_tests
    PRIVATE
        mpmc::mesh_petsc
        mpmc::discretization_petsc
        mpmc::flow_discretization)
set_target_properties(
    mpmc_mesh_discretization_petsc_tests
    PROPERTIES CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_mesh_discretization_petsc_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

add_test(
    NAME mesh.petsc.synthetic_2rank
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:mpmc_mesh_discretization_petsc_tests>
        ${MPIEXEC_POSTFLAGS})
set_tests_properties(
    mesh.petsc.synthetic_2rank
    PROPERTIES LABELS "mesh;mesh.petsc;integration;mpi" TIMEOUT 60)

add_test(
    NAME discretization.petsc.sparsity_owner_targeted_3rank
    COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 3
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:mpmc_mesh_discretization_petsc_tests>
        ${MPIEXEC_POSTFLAGS})
set_tests_properties(
    discretization.petsc.sparsity_owner_targeted_3rank
    PROPERTIES
        LABELS "discretization;discretization.petsc;integration;mpi;sparsity"
        TIMEOUT 60
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
