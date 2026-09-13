cmake_minimum_required(VERSION 3.24)
project(mpmc_pt_installer LANGUAGES NONE)

set(MPMC_PT_STAGING_DIR "" CACHE PATH
    "Absolute path to one validated PT product staging tree")
set(MPMC_PT_INSTALLER_TARGET "" CACHE STRING
    "Installer target: linux-x86_64, windows-x86_64, or macos-arm64")

if(NOT IS_ABSOLUTE "${MPMC_PT_STAGING_DIR}" OR
   NOT IS_DIRECTORY "${MPMC_PT_STAGING_DIR}")
    message(FATAL_ERROR
        "MPMC_PT_STAGING_DIR must name an existing absolute staging tree")
endif()

set(mpmc_staging_manifest
    "${MPMC_PT_STAGING_DIR}/share/mpmc-pt/product-staging-manifest.json")
if(NOT EXISTS "${mpmc_staging_manifest}")
    message(FATAL_ERROR "The PT product staging manifest is missing")
endif()
file(READ "${mpmc_staging_manifest}" mpmc_staging_json)
string(JSON mpmc_staging_convention GET "${mpmc_staging_json}" convention)
string(JSON mpmc_product_id GET "${mpmc_staging_json}" product id)
string(JSON mpmc_product_version GET "${mpmc_staging_json}" product version)
string(JSON mpmc_dependency_profile
    GET "${mpmc_staging_json}" platform dependency_profile)
string(JSON mpmc_entry_path GET "${mpmc_staging_json}" entry path)

if(NOT mpmc_staging_convention STREQUAL "MPMC/PT/product-staging/v1" OR
   NOT mpmc_product_id STREQUAL "mpmc-pt")
    message(FATAL_ERROR "The installer input is not a PT product staging v1 tree")
endif()
if(NOT mpmc_product_version MATCHES
   "^[0-9]+[.][0-9]+[.][0-9]+([-+][A-Za-z0-9.-]+)?$")
    message(FATAL_ERROR "The staged PT product version is not package-safe")
endif()
if(IS_ABSOLUTE "${mpmc_entry_path}" OR
   mpmc_entry_path MATCHES "(^|/)\\.\\.(/|$)" OR
   NOT EXISTS "${MPMC_PT_STAGING_DIR}/${mpmc_entry_path}")
    message(FATAL_ERROR "The staged PT entry path is invalid")
endif()

if(MPMC_PT_INSTALLER_TARGET STREQUAL "linux-x86_64")
    set(mpmc_expected_profile "linux-x86_64-release")
    set(CPACK_GENERATOR "DEB")
    set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/mpmc-pt")
    set(CPACK_DEBIAN_PACKAGE_NAME "mpmc-pt")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "MPMC_HNU")
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
    set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS OFF)
elseif(MPMC_PT_INSTALLER_TARGET STREQUAL "windows-x86_64")
    set(mpmc_expected_profile "windows-x86_64-release")
    set(CPACK_GENERATOR "WIX")
    set(CPACK_WIX_VERSION 3)
    set(CPACK_WIX_ARCHITECTURE "x64")
    set(CPACK_WIX_INSTALL_SCOPE "perMachine")
    set(CPACK_WIX_ROOT_FOLDER_ID "ProgramFiles64Folder")
    # UUIDv5(NAMESPACE_URL,
    # "https://github.com/lihuihnu/MPMC_HNU/products/pt") is stable for
    # upgrades and belongs only to this product identity.
    set(CPACK_WIX_UPGRADE_GUID "4D54F0D2-922A-5E2C-98E4-F0640F891BBA")
    # Project licensing is not yet decided.  Display the explicit candidate
    # notice instead of CMake's generic license text in the WiX UI.
    set(CPACK_RESOURCE_FILE_LICENSE
        "${CMAKE_CURRENT_LIST_DIR}/candidate-notice.txt")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "MPMC PT")
elseif(MPMC_PT_INSTALLER_TARGET STREQUAL "macos-arm64")
    set(mpmc_expected_profile "macos-armv8-release")
    set(CPACK_GENERATOR "productbuild")
    set(CPACK_PACKAGING_INSTALL_PREFIX
        "/Library/Application Support/MPMC/PT")
    set(CPACK_PRODUCTBUILD_IDENTIFIER "io.github.lihuihnu.mpmc.pt")
    set(CPACK_PRODUCTBUILD_DOMAINS TRUE)
    set(CPACK_PRODUCTBUILD_DOMAINS_ANYWHERE FALSE)
    set(CPACK_PRODUCTBUILD_DOMAINS_USER FALSE)
    set(CPACK_PRODUCTBUILD_DOMAINS_ROOT TRUE)
else()
    message(FATAL_ERROR
        "MPMC_PT_INSTALLER_TARGET must be linux-x86_64, "
        "windows-x86_64, or macos-arm64")
endif()

if(NOT mpmc_dependency_profile STREQUAL mpmc_expected_profile)
    message(FATAL_ERROR
        "The staging dependency profile does not match the installer target")
endif()

# The installer consumes the completed staging tree as opaque payload.  It
# never links the PT service or reaches into any model-specific implementation.
install(DIRECTORY "${MPMC_PT_STAGING_DIR}/"
    DESTINATION "."
    USE_SOURCE_PERMISSIONS)
install(FILES "${CMAKE_CURRENT_LIST_DIR}/candidate-notice.txt"
    DESTINATION "share/mpmc-pt")

set(CPACK_PACKAGE_NAME "MPMC PT Runtime")
set(CPACK_PACKAGE_VENDOR "MPMC_HNU")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Model-neutral native PT service runtime")
set(CPACK_PACKAGE_VERSION "${mpmc_product_version}")
set(CPACK_PACKAGE_FILE_NAME
    "mpmc-pt-${mpmc_product_version}-${MPMC_PT_INSTALLER_TARGET}")
set(CPACK_PACKAGE_CHECKSUM "SHA256")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY FALSE)
set(CPACK_VERBATIM_VARIABLES TRUE)

include(CPack)
