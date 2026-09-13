cmake_minimum_required(VERSION 3.21)
project(mpmc_pt_product VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
include(GNUInstallDirs)

set(MPMC_PT_PRODUCT_BUILD_REVISION "unversioned-local-build" CACHE STRING
    "Stable source revision recorded in the staged product manifest")
set(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER "" CACHE STRING
    "Dependency provider used for this product staging build")

if(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER STREQUAL "" AND
   DEFINED MPMC_PT_PRODUCT_DEPENDENCY_METADATA_DIR)
    set(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER "conan")
elseif(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER STREQUAL "" AND
       DEFINED VCPKG_TARGET_TRIPLET)
    set(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER "vcpkg")
endif()

if(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER STREQUAL "conan")
    foreach(required_variable IN ITEMS
            MPMC_PT_PRODUCT_DEPENDENCY_METADATA_DIR
            MPMC_PT_PRODUCT_OPENSSL_EXECUTABLE
            MPMC_PT_PRODUCT_CONAN_PROFILE_ID
            MPMC_PT_PRODUCT_CONAN_LOCK_SHA256
            MPMC_PT_PRODUCT_CONAN_VERSION)
        if(NOT DEFINED ${required_variable} OR
           "${${required_variable}}" STREQUAL "")
            message(FATAL_ERROR
                "The PT product staging build requires ${required_variable} "
                "from a validated restore-only Conan graph")
        endif()
    endforeach()
    if(NOT MPMC_PT_PRODUCT_CONAN_PROFILE_ID MATCHES
           "^[A-Za-z0-9][A-Za-z0-9._+-]*$")
        message(FATAL_ERROR "The Conan profile id is not a stable identifier")
    endif()
    string(LENGTH "${MPMC_PT_PRODUCT_CONAN_LOCK_SHA256}"
           mpmc_conan_lock_sha256_length)
    if(NOT mpmc_conan_lock_sha256_length EQUAL 64 OR
       NOT MPMC_PT_PRODUCT_CONAN_LOCK_SHA256 MATCHES "^[0-9a-f]+$")
        message(FATAL_ERROR "The Conan lock digest must be a SHA-256 value")
    endif()
    if(NOT EXISTS
       "${MPMC_PT_PRODUCT_DEPENDENCY_METADATA_DIR}/dependency-manifest.json")
        message(FATAL_ERROR "The generated Conan dependency manifest is missing")
    endif()
    if(NOT IS_DIRECTORY
       "${MPMC_PT_PRODUCT_DEPENDENCY_METADATA_DIR}/third-party")
        message(FATAL_ERROR "The generated Conan license set is missing")
    endif()
    if(NOT EXISTS "${MPMC_PT_PRODUCT_OPENSSL_EXECUTABLE}")
        message(FATAL_ERROR "The OpenSSL executable for the mTLS smoke is missing")
    endif()
    set(MPMC_PT_PRODUCT_DEPENDENCY_PROFILE
        "${MPMC_PT_PRODUCT_CONAN_PROFILE_ID}")
    set(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER_VERSION
        "${MPMC_PT_PRODUCT_CONAN_VERSION}")
    set(MPMC_PT_PRODUCT_DEPENDENCY_LOCK_KIND "conan-lock-sha256")
    set(MPMC_PT_PRODUCT_DEPENDENCY_LOCK_ID
        "${MPMC_PT_PRODUCT_CONAN_LOCK_SHA256}")
elseif(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER STREQUAL "vcpkg")
    if(NOT DEFINED VCPKG_TARGET_TRIPLET OR
       NOT DEFINED VCPKG_INSTALLED_DIR)
        message(FATAL_ERROR
            "The vcpkg target triplet and installed directory are required")
    endif()
    file(READ "${CMAKE_CURRENT_LIST_DIR}/vcpkg.json" mpmc_vcpkg_manifest)
    string(JSON MPMC_PT_VCPKG_BASELINE
        GET "${mpmc_vcpkg_manifest}" builtin-baseline)
    string(LENGTH "${MPMC_PT_VCPKG_BASELINE}" mpmc_vcpkg_baseline_length)
    if(NOT mpmc_vcpkg_baseline_length EQUAL 40 OR
       NOT MPMC_PT_VCPKG_BASELINE MATCHES "^[0-9a-f]+$")
        message(FATAL_ERROR "The PT product vcpkg baseline must be a full SHA")
    endif()
    set(MPMC_PT_PRODUCT_DEPENDENCY_PROFILE "${VCPKG_TARGET_TRIPLET}")
    set(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER_VERSION "builtin")
    set(MPMC_PT_PRODUCT_DEPENDENCY_LOCK_KIND "vcpkg-builtin-baseline")
    set(MPMC_PT_PRODUCT_DEPENDENCY_LOCK_ID "${MPMC_PT_VCPKG_BASELINE}")
else()
    message(FATAL_ERROR
        "MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER must be conan or vcpkg")
endif()
string(LENGTH "${MPMC_PT_PRODUCT_BUILD_REVISION}"
       mpmc_product_build_revision_length)
if(mpmc_product_build_revision_length LESS 1 OR
   mpmc_product_build_revision_length GREATER 128 OR
   NOT MPMC_PT_PRODUCT_BUILD_REVISION MATCHES
       "^[A-Za-z0-9][A-Za-z0-9._+-]*$")
    message(FATAL_ERROR
        "MPMC_PT_PRODUCT_BUILD_REVISION must be a stable non-secret identifier")
endif()

get_filename_component(mpmc_repository_root
    "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
add_subdirectory("${mpmc_repository_root}/modules/pt_process"
                 pt-process-library)

function(mpmc_pt_product_enable_strict_warnings target)
    if(MSVC)
        # Imported dependencies are SYSTEM by default; INCLUDE holds the MSVC
        # and Windows SDK roots. Keep project /I paths at /W4 /WX while
        # preventing external-header diagnostics from being promoted by /WX.
        target_compile_options(${target} PRIVATE
            /W4 /WX /permissive- /utf-8
            /external:env:INCLUDE /external:W0)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror)
    endif()
endfunction()

foreach(target IN ITEMS
        mpmc_runtime_grpc
        mpmc_pt_process
        mpmc_pt_service_host)
    mpmc_pt_product_enable_strict_warnings(${target})
endforeach()

if(APPLE)
    set_target_properties(mpmc_pt_service_host PROPERTIES
        INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
    set_target_properties(mpmc_pt_service_host PROPERTIES
        INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

set(mpmc_runtime_dependency_filters
    PRE_EXCLUDE_REGEXES
        "^api-ms-win-.*"
        "^ext-ms-win-.*")
if(WIN32)
    list(APPEND mpmc_runtime_dependency_filters
        POST_EXCLUDE_REGEXES
            ".*[/\\\\]Windows[/\\\\]System32[/\\\\].*")
elseif(APPLE)
    list(APPEND mpmc_runtime_dependency_filters
        POST_EXCLUDE_REGEXES
            "^/System/Library/.*"
            "^/usr/lib/.*")
elseif(UNIX)
    list(APPEND mpmc_runtime_dependency_filters
        POST_EXCLUDE_REGEXES
            "^/lib/.*"
            "^/lib64/.*"
            "^/usr/lib/.*"
            "^/usr/lib64/.*")
endif()

install(TARGETS mpmc_pt_service_host
    RUNTIME_DEPENDENCIES ${mpmc_runtime_dependency_filters}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    FRAMEWORK DESTINATION ".")

configure_file(
    "${CMAKE_CURRENT_LIST_DIR}/product-staging-manifest.json.in"
    "${CMAKE_CURRENT_BINARY_DIR}/product-staging-manifest.json"
    @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/product-staging-manifest.json"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/mpmc-pt")
if(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER STREQUAL "conan")
    install(FILES
        "${MPMC_PT_PRODUCT_DEPENDENCY_METADATA_DIR}/dependency-manifest.json"
        DESTINATION "${CMAKE_INSTALL_DATADIR}/mpmc-pt"
        RENAME "dependency-manifest.json")
    install(DIRECTORY
        "${MPMC_PT_PRODUCT_DEPENDENCY_METADATA_DIR}/third-party/"
        DESTINATION "${CMAKE_INSTALL_DATADIR}/mpmc-pt/third-party")
else()
    install(FILES "${CMAKE_CURRENT_LIST_DIR}/vcpkg.json"
        DESTINATION "${CMAKE_INSTALL_DATADIR}/mpmc-pt"
        RENAME "dependency-manifest.json")
    set(mpmc_vcpkg_share
        "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
    foreach(package IN ITEMS grpc protobuf openssl)
        if(NOT EXISTS "${mpmc_vcpkg_share}/${package}/copyright")
            message(FATAL_ERROR
                "Locked vcpkg dependency ${package} has no copyright file")
        endif()
    endforeach()
    file(GLOB mpmc_dependency_copyrights CONFIGURE_DEPENDS
        "${mpmc_vcpkg_share}/*/copyright")
    if(NOT mpmc_dependency_copyrights)
        message(FATAL_ERROR "The locked vcpkg dependency copyright set is empty")
    endif()
    foreach(copyright_file IN LISTS mpmc_dependency_copyrights)
        get_filename_component(package_directory "${copyright_file}" DIRECTORY)
        get_filename_component(package_name "${package_directory}" NAME)
        install(FILES "${copyright_file}"
            DESTINATION
                "${CMAKE_INSTALL_DATADIR}/mpmc-pt/third-party/${package_name}")
    endforeach()
endif()

if(BUILD_TESTING)
    find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
    if(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER STREQUAL "vcpkg")
        find_program(MPMC_PT_PRODUCT_OPENSSL_EXECUTABLE NAMES openssl
            HINTS
                "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/tools/openssl"
            REQUIRED)
    endif()
    add_executable(mpmc_pt_product_probe
        "${CMAKE_CURRENT_LIST_DIR}/tests/staged_probe.cpp")
    target_compile_features(mpmc_pt_product_probe PRIVATE cxx_std_20)
    target_link_libraries(mpmc_pt_product_probe PRIVATE mpmc::runtime_grpc)
    set_target_properties(mpmc_pt_product_probe PROPERTIES CXX_EXTENSIONS OFF)
    mpmc_pt_product_enable_strict_warnings(mpmc_pt_product_probe)

    add_executable(mpmc_pt_product_contract_tests
        "${mpmc_repository_root}/tests/runtime/pt_grpc_adapter/pt_grpc_adapter_test.cpp"
        "${mpmc_repository_root}/tests/runtime/pt_grpc_adapter/pt_grpc_adapter_headers.cpp")
    target_compile_features(mpmc_pt_product_contract_tests PRIVATE cxx_std_20)
    target_link_libraries(mpmc_pt_product_contract_tests
        PRIVATE mpmc::runtime_grpc)
    set_target_properties(mpmc_pt_product_contract_tests
        PROPERTIES CXX_EXTENSIONS OFF)
    mpmc_pt_product_enable_strict_warnings(mpmc_pt_product_contract_tests)

    foreach(case IN ITEMS
            mapping_and_outcomes
            wire_and_size_limits
            concurrency_gate
            deadline_and_cancel
            headers)
        add_test(NAME product.pt.contract.${case}
            COMMAND mpmc_pt_product_contract_tests ${case})
        set_tests_properties(product.pt.contract.${case} PROPERTIES
            LABELS "product;staging;runtime;grpc;contract"
            TIMEOUT 30)
    endforeach()

    add_test(NAME product.pt.architecture_boundary
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_LIST_DIR}/tests/test_product_boundary.py")
    set_tests_properties(product.pt.architecture_boundary PROPERTIES
        LABELS "product;staging;architecture"
        TIMEOUT 30)

    add_test(NAME product.pt.restore_only_dependency_metadata
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_LIST_DIR}/tests/test_conan_metadata.py")
    set_tests_properties(product.pt.restore_only_dependency_metadata PROPERTIES
        LABELS "product;staging;dependencies;provenance"
        TIMEOUT 30)

    add_test(NAME product.pt.staged_mtls_discovery
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_LIST_DIR}/tests/staged_smoke.py"
            --cmake "${CMAKE_COMMAND}"
            --build-directory "${CMAKE_BINARY_DIR}"
            --stage-directory "${CMAKE_BINARY_DIR}/stage"
            --configuration "$<CONFIG>"
            --host-relative-path
                "${CMAKE_INSTALL_BINDIR}/mpmc_pt_service_host${CMAKE_EXECUTABLE_SUFFIX}"
            --probe "$<TARGET_FILE:mpmc_pt_product_probe>"
            --openssl "${MPMC_PT_PRODUCT_OPENSSL_EXECUTABLE}")
    set_tests_properties(product.pt.staged_mtls_discovery PROPERTIES
        LABELS "product;staging;runtime;grpc;mtls;provenance"
        TIMEOUT 60)
endif()
