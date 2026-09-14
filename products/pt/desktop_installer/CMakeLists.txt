cmake_minimum_required(VERSION 3.24)
project(mpmc_pt_windows_desktop_installer LANGUAGES NONE)

set(MPMC_PT_DESKTOP_PACKAGE_DIR "" CACHE PATH
    "Absolute path to the packaged Windows Electron application directory")
set(MPMC_PT_DESKTOP_EXPECTED_REVISION "" CACHE STRING
    "Exact source revision expected in the desktop preview manifest")

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

set(mpmc_desktop_executable
    "${MPMC_PT_DESKTOP_PACKAGE_DIR}/MPMC-PT-Desktop-Preview.exe")
set(mpmc_desktop_resources
    "${MPMC_PT_DESKTOP_PACKAGE_DIR}/resources")
set(mpmc_desktop_manifest
    "${mpmc_desktop_resources}/desktop-preview-manifest.json")
set(mpmc_desktop_native_manifest
    "${mpmc_desktop_resources}/desktop-native/share/mpmc-pt/product-staging-manifest.json")
set(mpmc_desktop_notice
    "${mpmc_desktop_resources}/desktop-preview-notice.txt")

foreach(required_path IN ITEMS
        "${mpmc_desktop_executable}"
        "${mpmc_desktop_manifest}"
        "${mpmc_desktop_native_manifest}"
        "${mpmc_desktop_notice}")
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

if(NOT mpmc_desktop_convention STREQUAL "MPMC/PT/desktop-preview/v1" OR
   NOT mpmc_desktop_revision STREQUAL "${MPMC_PT_DESKTOP_EXPECTED_REVISION}" OR
   NOT mpmc_desktop_platform STREQUAL "win32" OR
   NOT mpmc_desktop_architecture STREQUAL "x64" OR
   mpmc_desktop_release_eligible OR
   NOT mpmc_desktop_session_transport STREQUAL "native-grpc-loopback-bearer")
    message(FATAL_ERROR
        "Desktop preview manifest is incompatible with the Windows MSI gate")
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
set(CPACK_WIX_UPGRADE_GUID "40656E02-E6EE-59D5-A5A6-73F3C3C452D6")

set(CPACK_PACKAGE_NAME "MPMC PT Desktop")
set(CPACK_PACKAGE_VENDOR "MPMC_HNU")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "MPMC PT Electron desktop with embedded native PT runtime")
set(CPACK_PACKAGE_VERSION "0.1.0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "MPMC PT Desktop")
set(CPACK_PACKAGE_FILE_NAME
    "mpmc-pt-desktop-0.1.0-windows-x86_64")
set(CPACK_PACKAGE_CHECKSUM "SHA256")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY FALSE)
set(CPACK_VERBATIM_VARIABLES TRUE)
set(CPACK_RESOURCE_FILE_LICENSE
    "${CMAKE_CURRENT_LIST_DIR}/candidate-notice.txt")
set(CPACK_RESOURCE_FILE_README
    "${CMAKE_CURRENT_LIST_DIR}/candidate-notice.txt")
set(CPACK_RESOURCE_FILE_WELCOME
    "${CMAKE_CURRENT_LIST_DIR}/candidate-notice.txt")
set(CPACK_PACKAGE_EXECUTABLES
    "MPMC-PT-Desktop-Preview" "MPMC PT Desktop")
set(CPACK_CREATE_DESKTOP_LINKS "MPMC-PT-Desktop-Preview")

include(CPack)
