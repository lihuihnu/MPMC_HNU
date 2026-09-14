cmake_minimum_required(VERSION 3.24)
project(mpmc_pt_android_portability LANGUAGES CXX)

if(NOT ANDROID)
    message(FATAL_ERROR
        "products/pt_android is an Android-NDK-only portability target")
endif()

get_filename_component(mpmc_repository_root
    "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)

if(NOT TARGET mpmc::runtime)
    add_subdirectory("${mpmc_repository_root}/modules/runtime"
                     runtime-library)
endif()

# Reuse the canonical repository-curated parameter/backend assembly source
# directly. The Android-private composition_root.hpp shim exposes only the
# ownership primitive those translation units require, so runtime_grpc and the
# process host never enter this target.
set(mpmc_android_canonical_pt_sources
    "${mpmc_repository_root}/modules/pt_process/src/configured_backends.cpp"
    "${mpmc_repository_root}/modules/pt_process/src/parameter_snapshot_supplier.cpp")

add_library(mpmc_pt_android_core SHARED
    src/portability.cpp
    src/jni_smoke.cpp
    src/jni_product.cpp
    ${mpmc_android_canonical_pt_sources})

target_compile_features(mpmc_pt_android_core PRIVATE cxx_std_20)
target_include_directories(mpmc_pt_android_core
    PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/include"
    PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/private_include"
        "${mpmc_repository_root}/modules/pt_process/include")
target_link_libraries(mpmc_pt_android_core PRIVATE mpmc::runtime)
set_target_properties(mpmc_pt_android_core PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
    OUTPUT_NAME "mpmc_pt_android_core")

target_compile_options(mpmc_pt_android_core PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
