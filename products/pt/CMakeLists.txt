cmake_minimum_required(VERSION 3.21)
project(mpmc_pt_product VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
include(GNUInstallDirs)

option(MPMC_PT_PRODUCT_REQUIRE_VCPKG
       "Require the locked product vcpkg manifest" ON)
set(MPMC_PT_PRODUCT_BUILD_REVISION "unversioned-local-build" CACHE STRING
    "Stable source revision recorded in the staged product manifest")

if(MPMC_PT_PRODUCT_REQUIRE_VCPKG AND NOT DEFINED VCPKG_TARGET_TRIPLET)
    message(FATAL_ERROR
        "The PT product staging build requires the locked vcpkg manifest. "
        "Configure with "
        "-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake")
endif()

file(READ "${CMAKE_CURRENT_LIST_DIR}/vcpkg.json" mpmc_vcpkg_manifest)
string(JSON MPMC_PT_VCPKG_BASELINE
    GET "${mpmc_vcpkg_manifest}" builtin-baseline)
string(LENGTH "${MPMC_PT_VCPKG_BASELINE}" mpmc_vcpkg_baseline_length)
if(NOT mpmc_vcpkg_baseline_length EQUAL 40 OR
   NOT MPMC_PT_VCPKG_BASELINE MATCHES "^[0-9a-f]+$")
    message(FATAL_ERROR "The PT product vcpkg baseline must be a full commit SHA")
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
        target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
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
install(FILES "${CMAKE_CURRENT_LIST_DIR}/vcpkg.json"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/mpmc-pt"
    RENAME "dependency-manifest.json")

if(DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
    set(mpmc_vcpkg_share
        "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/share")
    foreach(package IN ITEMS grpc protobuf openssl)
        if(NOT EXISTS "${mpmc_vcpkg_share}/${package}/copyright")
            message(FATAL_ERROR
                "Locked dependency ${package} has no staged copyright file")
        endif()
    endforeach()
    file(GLOB mpmc_dependency_copyrights CONFIGURE_DEPENDS
        "${mpmc_vcpkg_share}/*/copyright")
    if(NOT mpmc_dependency_copyrights)
        message(FATAL_ERROR "The locked dependency copyright set is empty")
    endif()
    foreach(copyright_file IN LISTS mpmc_dependency_copyrights)
        get_filename_component(package_directory "${copyright_file}" DIRECTORY)
        get_filename_component(package_name "${package_directory}" NAME)
        install(FILES "${copyright_file}"
            DESTINATION
                "${CMAKE_INSTALL_DATADIR}/mpmc-pt/third-party/${package_name}")
    endforeach()
elseif(MPMC_PT_PRODUCT_REQUIRE_VCPKG)
    message(FATAL_ERROR "The locked vcpkg installation directory is unavailable")
endif()

if(BUILD_TESTING)
    find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
    find_program(MPMC_PT_PRODUCT_OPENSSL_EXECUTABLE NAMES openssl
        HINTS
            "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/tools/openssl"
        REQUIRED)

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
