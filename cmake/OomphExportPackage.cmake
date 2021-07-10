# =============================================================================
# Create config files that can be included by other projects to find and use
# oomph-lib using find_package(...) or find_package(... CONFIG).
# =============================================================================
# Introduce the variables CMAKE_INSTALL_{LIBDIR,BINDIR,INCLUDEDIR} and functions
# write_basic_package_version_file(...) and configure_package_config_file(...)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ------------------------------------------------------------------------------
# Define additional variables required by the configure file here

# Get the compile definitions
get_directory_property(
  OOMPH_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                      COMPILE_DEFINITIONS)

# ------------------------------------------------------------------------------

# Configure 'oomphlibConfig.cmake' using 'oomphlibConfig.cmake.in'
configure_package_config_file(
  ${OOMPH_CMAKE_CONFIG_TEMPLATE} ${OOMPH_CMAKE_PROJECT_CONFIG_FILE}
  INSTALL_DESTINATION ${OOMPH_CONFIG_INSTALL_DIR})

# Configure 'oomphlibConfigVersion.cmake' which exports the version info
write_basic_package_version_file(
  ${OOMPH_CMAKE_VERSION_CONFIG_FILE}
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY SameMajorVersion)

# Add all targets to the build-tree export set
export(
  EXPORT ${OOMPH_TARGETS_EXPORT_NAME}
  NAMESPACE ${PROJECT_NAMESPACE}::
  FILE ${OOMPH_CMAKE_PROJECT_EXPORTS_FILE})

# Install config files
install(FILES ${OOMPH_CMAKE_PROJECT_CONFIG_FILE}
              ${OOMPH_CMAKE_VERSION_CONFIG_FILE}
        DESTINATION ${OOMPH_CONFIG_INSTALL_DIR})

# Install the export targets file
install(
  EXPORT ${OOMPH_TARGETS_EXPORT_NAME}
  NAMESPACE ${PROJECT_NAMESPACE}::
  FILE "${OOMPH_TARGETS_EXPORT_NAME}.cmake"
  DESTINATION ${OOMPH_CONFIG_INSTALL_DIR})

# Install pkg-config
configure_file("${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}.pc.in"
               "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Handle everything else that should be exported, e.g. useful CMake modules:

# The list of modules to export
set(CMAKE_MODULES_TO_INSTALL OomphAddExecutable.cmake OomphAddTest.cmake)

# Place each module in the appropriate build/install directory
foreach(MODULE IN LISTS CMAKE_MODULES_TO_INSTALL)
  # Copy to binary (i.e. build) directory
  file(COPY "${CMAKE_CURRENT_LIST_DIR}/${MODULE}"
       DESTINATION ${OOMPH_CMAKE_CONFIG_DIR})

  # Install it
  install(FILES "${CMAKE_CURRENT_LIST_DIR}/${MODULE}"
          DESTINATION ${OOMPH_CONFIG_INSTALL_DIR})
endforeach()
# ------------------------------------------------------------------------------
