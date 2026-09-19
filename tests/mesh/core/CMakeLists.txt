cmake_minimum_required(VERSION 3.21)
project(mpmc_mesh_core_tests LANGUAGES CXX)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../../modules/mesh" mesh-library)
enable_testing()

add_executable(mpmc_mesh_core_tests
    core_test.cpp
    entity_header_self_contained.cpp
    csr_header_self_contained.cpp
    corner_point_geometry_3d_header_self_contained.cpp
    topology_header_self_contained.cpp
    cartesian_2d_header_self_contained.cpp
    geometry_2d_header_self_contained.cpp
    gmsh_4_1_header_self_contained.cpp
    grdecl_header_self_contained.cpp
    vtu_header_self_contained.cpp
    face_boundary_header_self_contained.cpp
    dense_field_header_self_contained.cpp
    dof_layout_header_self_contained.cpp
    partition_snapshot_header_self_contained.cpp
    dof_numbering_header_self_contained.cpp
    shared_entity_plan_header_self_contained.cpp)
target_link_libraries(mpmc_mesh_core_tests PRIVATE mpmc::mesh)
set_target_properties(mpmc_mesh_core_tests PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(mpmc_mesh_core_tests PRIVATE /W4 /WX /permissive- /utf-8)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(mpmc_mesh_core_tests PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
endif()

option(MPMC_MESH_ENABLE_SANITIZERS
       "Use AddressSanitizer and UndefinedBehaviorSanitizer for mesh core tests" OFF)
if(MPMC_MESH_ENABLE_SANITIZERS)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR
       NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        message(FATAL_ERROR "Mesh core sanitizers require Linux GCC/Clang")
    endif()
    target_compile_options(mpmc_mesh_core_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
    target_link_options(mpmc_mesh_core_tests PRIVATE
        -fsanitize=address,undefined -fno-sanitize-recover=all)
endif()

foreach(case IN ITEMS strong_indices two_cell_quad invalid_csr topology_snapshot topology_invalid cartesian_2d_topology cartesian_2d_invalid cartesian_2d_geometry cartesian_2d_geometry_invalid face_boundary_snapshot face_boundary_invalid dense_field_snapshot dense_field_invalid gmsh_4_1_import gmsh_4_1_roundtrip gmsh_4_1_invalid vtu_ascii_import vtu_ascii_roundtrip vtu_ascii_invalid grdecl_import grdecl_inactive_degenerate grdecl_invalid dof_layout_snapshot dof_layout_invalid partition_serial_snapshot partition_local_snapshot partition_invalid dof_numbering_serial dof_numbering_local dof_numbering_invalid shared_entity_plan_pair shared_entity_plan_invalid headers)
    add_test(NAME mesh.core.${case} COMMAND mpmc_mesh_core_tests ${case})
    set_tests_properties(mesh.core.${case}
        PROPERTIES LABELS "mesh;mesh.core;unit" TIMEOUT 30)
endforeach()
