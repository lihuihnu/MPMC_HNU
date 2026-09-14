cmake_minimum_required(VERSION 3.24)
project(mpmc_pt_windows_desktop_installer LANGUAGES NONE)

set(MPMC_PT_DESKTOP_PACKAGE_DIR "" CACHE PATH
    "Absolute path to the packaged Windows Electron application directory")
set(MPMC_PT_DESKTOP_EXPECTED_REVISION "" CACHE STRING
    "Exact source revision expected in the desktop package manifest")
set(MPMC_PT_DESKTOP_RELEASE_IDENTITY
    "${CMAKE_CURRENT_LIST_DIR}/release-identity.json" CACHE FILEPATH
    "Authoritative Windows desktop release identity contract")
option(MPMC_PT_DESKTOP_RELEASE_CANDIDATE
    "Build the fixed signed-release-candidate MSI identity" OFF)

if(NOT IS_ABSOLUTE "${MPMC_PT_DESKTOP_PACKAGE_DIR}" OR
   NOT IS_DIRECTORY "${MPMC_PT_DESKTOP_PACKAGE_DIR}")
    message(FATAL_ERROR
        "MPMC_PT_DESKTOP_PACKAGE_DIR must name an existing absolute directory")
endif()
string(LENGTH "${MPMC_PT_DESKTOP_EXPECTED_REVISION}" mpmc_revision_length)
if(NOT mpmc_revision_length EQUAL 40 OR
   NOT MPMC_PT_DESKTOP_EXPECTED_REVISION MATCHES "^[0-9a-f]+$")
    message(FATAL_ERROR
        "MPMC_PT_DESKTOP_EXPECTED_REVISION must be one full lowercase SHA")
endif()
if(NOT IS_ABSOLUTE "${MPMC_PT_DESKTOP_RELEASE_IDENTITY}" OR
   NOT EXISTS "${MPMC_PT_DESKTOP_RELEASE_IDENTITY}")
    message(FATAL_ERROR
        "MPMC_PT_DESKTOP_RELEASE_IDENTITY must name an existing absolute JSON file")
endif()

file(READ "${MPMC_PT_DESKTOP_RELEASE_IDENTITY}" mpmc_release_identity_json)
string(JSON mpmc_identity_convention GET "${mpmc_release_identity_json}" convention)
string(JSON mpmc_product_id GET "${mpmc_release_identity_json}" product_id)
string(JSON mpmc_product_name GET "${mpmc_release_identity_json}" product_name)
string(JSON mpmc_product_vendor GET "${mpmc_release_identity_json}" publisher)
string(JSON mpmc_identity_platform GET "${mpmc_release_identity_json}" platform)
string(JSON mpmc_install_directory GET "${mpmc_release_identity_json}" install_directory)
string(JSON mpmc_upgrade_guid GET "${mpmc_release_identity_json}" upgrade_guid)
string(JSON mpmc_release_channel GET "${mpmc_release_identity_json}" release channel)
string(JSON mpmc_release_sequence GET "${mpmc_release_identity_json}" release sequence)
string(JSON mpmc_display_version GET "${mpmc_release_identity_json}" release display_version)
string(JSON mpmc_msi_product_version GET "${mpmc_release_identity_json}" release msi_product_version)
string(JSON mpmc_rc_product_guid GET "${mpmc_release_identity_json}" release product_guid)
string(JSON mpmc_engineering_executable GET "${mpmc_release_identity_json}" application engineering_executable)
string(JSON mpmc_release_executable GET "${mpmc_release_identity_json}" application release_executable)
string(JSON mpmc_file_digest GET "${mpmc_release_identity_json}" authenticode file_digest)
string(JSON mpmc_timestamp_protocol GET "${mpmc_release_identity_json}" authenticode timestamp_protocol)
string(JSON mpmc_timestamp_digest GET "${mpmc_release_identity_json}" authenticode timestamp_digest)

set(mpmc_guid_regex "^[0-9A-F]{8}-[0-9A-F]{4}-[0-9A-F]{4}-[0-9A-F]{4}-[0-9A-F]{12}$")
if(NOT mpmc_identity_convention STREQUAL
       "MPMC/PT/windows-desktop-release-identity/v1" OR
   NOT mpmc_product_id STREQUAL "mpmc-pt-desktop" OR
   NOT mpmc_identity_platform STREQUAL "windows-x86_64" OR
   NOT mpmc_release_channel STREQUAL "rc" OR
   NOT mpmc_release_sequence EQUAL 1 OR
   NOT mpmc_display_version MATCHES "^[0-9]+[.][0-9]+[.][0-9]+-rc[.][1-9][0-9]*$" OR
   NOT mpmc_msi_product_version MATCHES "^[0-9]+[.][0-9]+[.][0-9]+$" OR
   NOT mpmc_upgrade_guid MATCHES "${mpmc_guid_regex}" OR
   NOT mpmc_rc_product_guid MATCHES "${mpmc_guid_regex}" OR
   NOT mpmc_file_digest STREQUAL "SHA256" OR
   NOT mpmc_timestamp_protocol STREQUAL "RFC3161" OR
   NOT mpmc_timestamp_digest STREQUAL "SHA256")
    message(FATAL_ERROR "The fixed Windows desktop release identity is invalid")
endif()

if(MPMC_PT_DESKTOP_RELEASE_CANDIDATE)
    set(mpmc_expected_desktop_convention
        "MPMC/PT/desktop-release-candidate-payload/v1")
    set(mpmc_desktop_executable_name "${mpmc_release_executable}")
    set(mpmc_desktop_notice_name "desktop-release-candidate-notice.txt")
    set(mpmc_installer_notice
        "${CMAKE_CURRENT_LIST_DIR}/release-candidate-notice.txt")
    set(mpmc_package_file_version "${mpmc_display_version}")
else()
    set(mpmc_expected_desktop_convention "MPMC/PT/desktop-preview/v1")
    set(mpmc_desktop_executable_name "${mpmc_engineering_executable}")
    set(mpmc_desktop_notice_name "desktop-preview-notice.txt")
    set(mpmc_installer_notice
        "${CMAKE_CURRENT_LIST_DIR}/candidate-notice.txt")
    set(mpmc_package_file_version "${mpmc_msi_product_version}")
endif()
string(REGEX REPLACE "[.]exe$" "" mpmc_desktop_executable_stem
    "${mpmc_desktop_executable_name}")

set(mpmc_desktop_executable
    "${MPMC_PT_DESKTOP_PACKAGE_DIR}/${mpmc_desktop_executable_name}")
set(mpmc_desktop_resources
    "${MPMC_PT_DESKTOP_PACKAGE_DIR}/resources")
set(mpmc_desktop_manifest
    "${mpmc_desktop_resources}/desktop-preview-manifest.json")
set(mpmc_desktop_native_manifest
    "${mpmc_desktop_resources}/desktop-native/share/mpmc-pt/product-staging-manifest.json")
set(mpmc_desktop_notice
    "${mpmc_desktop_resources}/${mpmc_desktop_notice_name}")

foreach(required_path IN ITEMS
        "${mpmc_desktop_executable}"
        "${mpmc_desktop_manifest}"
        "${mpmc_desktop_native_manifest}"
        "${mpmc_desktop_notice}"
        "${mpmc_installer_notice}")
    if(NOT EXISTS "${required_path}")
        message(FATAL_ERROR "Desktop installer input is incomplete: ${required_path}")
    endif()
endforeach()

file(READ "${mpmc_desktop_manifest}" mpmc_desktop_json)
string(JSON mpmc_desktop_convention GET "${mpmc_desktop_json}" convention)
string(JSON mpmc_desktop_revision GET "${mpmc_desktop_json}" source_revision)
string(JSON mpmc_desktop_platform GET "${mpmc_desktop_json}" platform)
string(JSON mpmc_desktop_architecture GET "${mpmc_desktop_json}" architecture)
string(JSON mpmc_desktop_release_eligible GET "${mpmc_desktop_json}" release_eligible)
string(JSON mpmc_desktop_session_transport
    GET "${mpmc_desktop_json}" desktop_session transport)

if(NOT mpmc_desktop_convention STREQUAL "${mpmc_expected_desktop_convention}" OR
   NOT mpmc_desktop_revision STREQUAL "${MPMC_PT_DESKTOP_EXPECTED_REVISION}" OR
   NOT mpmc_desktop_platform STREQUAL "win32" OR
   NOT mpmc_desktop_architecture STREQUAL "x64" OR
   mpmc_desktop_release_eligible OR
   NOT mpmc_desktop_session_transport STREQUAL "native-grpc-loopback-bearer")
    message(FATAL_ERROR
        "Desktop package manifest is incompatible with the Windows MSI gate")
endif()

if(MPMC_PT_DESKTOP_RELEASE_CANDIDATE)
    string(JSON mpmc_payload_is_rc GET "${mpmc_desktop_json}" release_candidate)
    string(JSON mpmc_payload_display_version GET "${mpmc_desktop_json}" release_identity display_version)
    string(JSON mpmc_payload_msi_version GET "${mpmc_desktop_json}" release_identity msi_product_version)
    string(JSON mpmc_payload_product_guid GET "${mpmc_desktop_json}" release_identity product_guid)
    string(JSON mpmc_payload_upgrade_guid GET "${mpmc_desktop_json}" release_identity upgrade_guid)
    if(NOT mpmc_payload_is_rc OR
       NOT mpmc_payload_display_version STREQUAL "${mpmc_display_version}" OR
       NOT mpmc_payload_msi_version STREQUAL "${mpmc_msi_product_version}" OR
       NOT mpmc_payload_product_guid STREQUAL "${mpmc_rc_product_guid}" OR
       NOT mpmc_payload_upgrade_guid STREQUAL "${mpmc_upgrade_guid}")
        message(FATAL_ERROR
            "Desktop RC payload does not match the fixed release identity")
    endif()
endif()

file(READ "${mpmc_desktop_native_manifest}" mpmc_native_json)
string(JSON mpmc_native_convention GET "${mpmc_native_json}" convention)
string(JSON mpmc_native_revision GET "${mpmc_native_json}" product build_revision)
string(JSON mpmc_native_profile GET "${mpmc_native_json}" platform dependency_profile)
string(JSON mpmc_native_transport GET "${mpmc_native_json}" entry transport)

if(NOT mpmc_native_convention STREQUAL "MPMC/PT/product-staging/v1" OR
   NOT mpmc_native_revision STREQUAL "${MPMC_PT_DESKTOP_EXPECTED_REVISION}" OR
   NOT mpmc_native_profile STREQUAL "windows-x86_64-release" OR
   NOT mpmc_native_transport STREQUAL "native-grpc-mtls")
    message(FATAL_ERROR
        "Embedded native PT staging does not match the desktop installer target")
endif()

install(DIRECTORY "${MPMC_PT_DESKTOP_PACKAGE_DIR}/"
    DESTINATION "."
    USE_SOURCE_PERMISSIONS)

set(CPACK_GENERATOR "WIX")
set(CPACK_WIX_VERSION 3)
set(CPACK_WIX_ARCHITECTURE "x64")
set(CPACK_WIX_INSTALL_SCOPE "perMachine")
set(CPACK_WIX_ROOT_FOLDER_ID "ProgramFiles64Folder")
set(CPACK_WIX_UPGRADE_GUID "${mpmc_upgrade_guid}")
if(MPMC_PT_DESKTOP_RELEASE_CANDIDATE)
    set(CPACK_WIX_PRODUCT_GUID "${mpmc_rc_product_guid}")
endif()

set(CPACK_PACKAGE_NAME "${mpmc_product_name}")
set(CPACK_PACKAGE_VENDOR "${mpmc_product_vendor}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "MPMC PT Electron desktop with embedded native PT runtime")
set(CPACK_PACKAGE_VERSION "${mpmc_msi_product_version}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${mpmc_install_directory}")
set(CPACK_PACKAGE_FILE_NAME
    "mpmc-pt-desktop-${mpmc_package_file_version}-windows-x86_64")
set(CPACK_PACKAGE_CHECKSUM "SHA256")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY FALSE)
set(CPACK_VERBATIM_VARIABLES TRUE)
set(CPACK_RESOURCE_FILE_LICENSE "${mpmc_installer_notice}")
set(CPACK_RESOURCE_FILE_README "${mpmc_installer_notice}")
set(CPACK_RESOURCE_FILE_WELCOME "${mpmc_installer_notice}")
set(CPACK_PACKAGE_EXECUTABLES
    "${mpmc_desktop_executable_stem}" "${mpmc_product_name}")
set(CPACK_CREATE_DESKTOP_LINKS "${mpmc_desktop_executable_stem}")

include(CPack)
