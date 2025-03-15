# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
# ...to be filled in...
#
# USAGE:
# ------
#
# ...to be filled in...
#
# EXAMPLE:
# --------
#
# ...to be filled in...
#
# =============================================================================
# cmake-format: on
include_guard()

# cmake-format: off
function(oomph_get_external_project_helper)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS CONFIGURE_HANDLED_BY_BUILD)
  set(SINGLE_VALUE_ARGS PROJECT_NAME URL INSTALL_DIR)
  set(MULTI_VALUE_ARGS
      PATCH_COMMAND
      CONFIGURE_COMMAND
      BUILD_COMMAND
      INSTALL_COMMAND
      TEST_COMMAND
      BUILD_BYPRODUCTS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the variables in this scope without a prefix for clarity
  set(PROJECT_NAME ${${PREFIX}_PROJECT_NAME})
  set(URL ${${PREFIX}_URL})
  set(INSTALL_DIR ${${PREFIX}_INSTALL_DIR})
  set(PATCH_COMMAND ${${PREFIX}_PATCH_COMMAND})
  set(CONFIGURE_COMMAND ${${PREFIX}_CONFIGURE_COMMAND})
  set(BUILD_COMMAND ${${PREFIX}_BUILD_COMMAND})
  set(INSTALL_COMMAND ${${PREFIX}_INSTALL_COMMAND})
  set(TEST_COMMAND ${${PREFIX}_TEST_COMMAND})
  set(BUILD_BYPRODUCTS ${${PREFIX}_BUILD_BYPRODUCTS})
  set(CONFIGURE_HANDLED_BY_BUILD FALSE)
  if(${${PREFIX}_CONFIGURE_HANDLED_BY_BUILD})
    set(CONFIGURE_HANDLED_BY_BUILD TRUE)
  endif()

  # If we want to disable tests, we just won't specify a test command to ExternalProject
  if(NOT OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS)
    set(TEST_COMMAND)
  endif()

  # Define how to configure/build/install the project
  ExternalProject_Add(
    ${PROJECT_NAME}
    INSTALL_DIR "${INSTALL_DIR}"
    URL "${URL}"
    LOG_DIR "${CMAKE_BINARY_DIR}/logs"
    BUILD_IN_SOURCE TRUE
    LOG_PATCH TRUE
    LOG_UPDATE TRUE
    LOG_DOWNLOAD TRUE
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
    LOG_TEST TRUE
    LOG_MERGED_STDOUTERR TRUE
    LOG_OUTPUT_ON_FAILURE TRUE
    UPDATE_DISCONNECTED TRUE
    BUILD_ALWAYS FALSE
    # Replace list separator ';' with a whitespace
    # LIST_SEPARATOR " "
    PATCH_COMMAND "${PATCH_COMMAND}"
    CONFIGURE_HANDLED_BY_BUILD ${CONFIGURE_HANDLED_BY_BUILD}
    CONFIGURE_COMMAND "${CONFIGURE_COMMAND}"
    BUILD_COMMAND "${BUILD_COMMAND}"
    TEST_COMMAND "${TEST_COMMAND}"
    INSTALL_COMMAND "${INSTALL_COMMAND}"
    BUILD_BYPRODUCTS "${BUILD_BYPRODUCTS}")
endfunction()
# cmake-format: on

# ---------------------------------------------------------------------------------
