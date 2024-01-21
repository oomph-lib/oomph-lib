# =============================================================================
# DESCRIPTION:
# ------------
# Create config files that can be included by other projects to find and use
# oomph-lib using find_package(...) or find_package(... CONFIG).
#
# NOTE: A number of the key OOMPH_... values referenced here will have been
# defined in OomphConfigureProject.cmake.
# =============================================================================
# Define additional variables required by the oomphlibConfig.cmake.in file here
include_guard()

# Get the compile definitions
get_directory_property(
  OOMPH_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                      COMPILE_DEFINITIONS)

# -----------------------------------------------------------------------------

# Introduce the functions write_basic_package_version_file(...) and
# configure_package_config_file(...)
include(CMakePackageConfigHelpers)

# Configure 'oomphlibConfig.cmake' using 'oomphlibConfig.cmake.in'
configure_package_config_file(
  "${OOMPH_CMAKE_CONFIG_FILE_TEMPLATE}" "${OOMPH_CMAKE_CONFIG_FILE}"
  INSTALL_DESTINATION "${OOMPH_INSTALL_CONFIG_DIR}"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Configure 'oomphlibConfigVersion.cmake' which exports the version info
write_basic_package_version_file(
  "${OOMPH_CMAKE_CONFIG_VERSION_FILE}"
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY SameMajorVersion)

# Add all targets to the build-tree export set
export(
  EXPORT ${OOMPH_EXPORTS_TARGET_NAME}
  NAMESPACE ${PROJECT_NAMESPACE}::
  FILE "${OOMPH_CMAKE_EXPORTS_FILE}")

# Install config files
install(FILES "${OOMPH_CMAKE_CONFIG_FILE}" "${OOMPH_CMAKE_CONFIG_VERSION_FILE}"
        DESTINATION "${OOMPH_INSTALL_CONFIG_DIR}")

# Install the export targets file
install(
  EXPORT ${OOMPH_EXPORTS_TARGET_NAME}
  NAMESPACE ${PROJECT_NAMESPACE}::
  FILE "${OOMPH_EXPORTS_TARGET_NAME}.cmake"
  DESTINATION "${OOMPH_INSTALL_CONFIG_DIR}")

# -----------------------------------------------------------------------------

# Install anything else that needs to be generated using configure_file(...),
# i.e. any .cmake.in files

# Install config files
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/oomphlibUninstall.cmake"
        DESTINATION "${OOMPH_INSTALL_CONFIG_DIR}")

# -----------------------------------------------------------------------------

# Handle all useful CMake modules that should be exported

# The list of modules to copy to the build directory
set(OOMPH_FILES_TO_COPY_TO_BUILD_DIR
    "${CMAKE_CURRENT_LIST_DIR}/OomphGetHashedTargetName.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphAddExecutable.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphAddTest.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphAddPureCppTest.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphDefineTestData.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphEnableCodeCoverage.cmake")

# The list of modules to distribute with the library
set(OOMPH_FILES_TO_INSTALL_TO_CMAKE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/FindGMP.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/FindMPFR.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/FindGKlib.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/FindMETIS.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/FindParMETIS.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/FindSuperLU_DIST.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphGetHashedTargetName.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphAddExecutable.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphAddTest.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphAddPureCppTest.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphDefineTestData.cmake"
    "${CMAKE_CURRENT_LIST_DIR}/OomphEnableCodeCoverage.cmake"
    "${OOMPH_ROOT_DIR}/scripts/fpdiff.py")

# Copy files to the build directory at configure-time
foreach(MODULE IN LISTS OOMPH_FILES_TO_COPY_TO_BUILD_DIR)
  file(COPY "${MODULE}" DESTINATION "${OOMPH_BUILD_DIR}")
endforeach()

# Copy files to the install directory at install-time
foreach(MODULE IN LISTS OOMPH_FILES_TO_INSTALL_TO_CMAKE_DIR)
  install(FILES "${MODULE}" DESTINATION "${OOMPH_INSTALL_CONFIG_DIR}")
endforeach()

# -----------------------------------------------------------------------------
