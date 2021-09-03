# ==============================================================================
# Modified version of the file here:
# https://github.com/pablospe/cmake-example-library/blob/master/cmake/SetEnv.cmake
# ==============================================================================
# Introduce CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_BINDIR, CMAKE_INSTALL_INCLUDEDIR
include(GNUInstallDirs)

# Important: prohibit in-source builds!
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# ------------------------------------------------------------------------------
# Configuration variables; create custom variables to avoid issues with regular
# variables being overriden by third-party libraries or by parties that import
# this library. The variable names used here follow a similar template to that
# used by the nlohmann JSON library:
# https://github.com/nlohmann/json/blob/develop/CMakeLists.txt
# ------------------------------------------------------------------------------

# The project directory; passed to consumers of the exported library incase they
# need to use anything from there, e.g. the demo_drivers nearly all need the
# fpdiff.py script from the scripts/ subdirectory.
set(OOMPH_ROOT_DIR "${PROJECT_SOURCE_DIR}")

# Set the install paths
set(OOMPH_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    CACHE INTERNAL "")
set(OOMPH_INSTALL_INCLUDE_DIR
    "${CMAKE_INSTALL_FULL_INCLUDEDIR}/${PROJECT_NAME}")
set(OOMPH_INSTALL_BIN_DIR "${CMAKE_INSTALL_FULL_BINDIR}/${PROJECT_NAME}")
set(OOMPH_INSTALL_LIB_DIR "${CMAKE_INSTALL_FULL_LIBDIR}/${PROJECT_NAME}")

# Define the export information
set(OOMPH_TARGETS_EXPORT_NAME "${PROJECT_NAME}Exports")
set(OOMPH_CMAKE_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(OOMPH_CMAKE_CONFIG_TEMPLATE
    "${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}Config.cmake.in")
set(OOMPH_CMAKE_VERSION_CONFIG_FILE
    "${OOMPH_CMAKE_CONFIG_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(OOMPH_CMAKE_PROJECT_CONFIG_FILE
    "${OOMPH_CMAKE_CONFIG_DIR}/${PROJECT_NAME}Config.cmake")
set(OOMPH_CMAKE_PROJECT_EXPORTS_FILE
    "${OOMPH_CMAKE_CONFIG_DIR}/${PROJECT_NAME}Exports.cmake")

# Silence warnings on MacOS about targets of ranlib having no symbols. This
# occurs when there is no code to be compiled, e.g. when #ifdef directives
# exclude all of the code from all of the files required by a library. This
# could possibly moved to it's own simple module. Found online somewhere...
if(APPLE)
  set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_C_ARCHIVE_FINISH
      "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  set(CMAKE_CXX_ARCHIVE_FINISH
      "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif(APPLE)

# Define the namespace for libraries to be exported within
set(PROJECT_NAMESPACE oomph)

# Storage for the list of libraries exported by oomph-lib
set(OOMPHLIB_LIBRARIES CACHE INTERNAL "" FORCE)

# Make sure different configurations (e.g. Debug and Release) don't collide
set(CMAKE_DEBUG_POSTFIX "d")

# Set the export target for libraries built by oomph-lib
set(TARGETS_EXPORT_NAME ${PROJECT_NAME}Exports)

# Define the 'uninstall' command to handle uninstalling the installed targets
configure_file(
  "${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Uninstall.cmake.in"
  "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake" IMMEDIATE @ONLY)
add_custom_target(
  uninstall COMMAND ${CMAKE_COMMAND} -P
                    ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake)
# ------------------------------------------------------------------------------
