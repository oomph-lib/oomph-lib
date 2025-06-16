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
  set(SINGLE_VALUE_ARGS PROJECT_NAME URL URL_HASH GIT_REPOSITORY GIT_TAG INSTALL_DIR GIT_SUBMODULES_RECURSE LIST_SEPARATOR)
  set(MULTI_VALUE_ARGS PATCH_COMMAND CMAKE_ARGS CONFIGURE_COMMAND BUILD_COMMAND INSTALL_COMMAND TEST_COMMAND BUILD_BYPRODUCTS GIT_SUBMODULES)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # If we want to disable tests, we just won't specify a test command to ExternalProject
  if(NOT OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS)
    set(ARG_TEST_COMMAND)
  endif()

  if(NOT ARG_GIT_REPOSITORY AND NOT ARG_URL)
    message(FATAL_ERROR "oomph_get_external_project: URL or GIT_REPOSITORY is required")
  endif()

  # Set default
  if(NOT ARG_GIT_SUBMODULES_RECURSE)
    set(ARG_GIT_SUBMODULES_RECURSE TRUE)
  endif()

  # Define how to configure/build/install the project
  if(ARG_GIT_REPOSITORY)
    ExternalProject_Add(
      ${ARG_PROJECT_NAME}
      GIT_REPOSITORY ${ARG_GIT_REPOSITORY}
      GIT_TAG ${ARG_GIT_TAG}
      GIT_SHALLOW TRUE
      GIT_SUBMODULES ${ARG_GIT_SUBMODULES}
      GIT_SUBMODULES_RECURSE ${ARG_GIT_SUBMODULES_RECURSE}
      LOG_DIR "${CMAKE_BINARY_DIR}/logs"
      INSTALL_DIR "${ARG_INSTALL_DIR}"
      BUILD_IN_SOURCE TRUE
      LIST_SEPARATOR ${ARG_LIST_SEPARATOR}
      CONFIGURE_COMMAND "${ARG_CONFIGURE_COMMAND}"
      BUILD_COMMAND "${ARG_BUILD_COMMAND}"
      INSTALL_COMMAND "${ARG_INSTALL_COMMAND}"
      PATCH_COMMAND "${ARG_PATCH_COMMAND}"
      TEST_COMMAND "${ARG_TEST_COMMAND}"
      CMAKE_ARGS "${ARG_CMAKE_ARGS}"
      BUILD_BYPRODUCTS "${ARG_BUILD_BYPRODUCTS}"
      CONFIGURE_HANDLED_BY_BUILD ${ARG_CONFIGURE_HANDLED_BY_BUILD}
      UPDATE_DISCONNECTED TRUE
      BUILD_ALWAYS FALSE
      LOG_DOWNLOAD ON
      LOG_UPDATE ON
      LOG_PATCH ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
      LOG_INSTALL ON
      LOG_TEST ON
      LOG_MERGED_STDOUTERR ON
      LOG_OUTPUT_ON_FAILURE ON)
  else()
    ExternalProject_Add(
      ${ARG_PROJECT_NAME}
      URL ${ARG_URL}
      URL_HASH ${ARG_URL_HASH}
      INSTALL_DIR "${ARG_INSTALL_DIR}"
      LOG_DIR "${CMAKE_BINARY_DIR}/logs"
      BUILD_IN_SOURCE TRUE
      LIST_SEPARATOR ${ARG_LIST_SEPARATOR}
      CONFIGURE_COMMAND "${ARG_CONFIGURE_COMMAND}"
      BUILD_COMMAND "${ARG_BUILD_COMMAND}"
      INSTALL_COMMAND "${ARG_INSTALL_COMMAND}"
      PATCH_COMMAND "${ARG_PATCH_COMMAND}"
      TEST_COMMAND "${ARG_TEST_COMMAND}"
      CMAKE_ARGS "${ARG_CMAKE_ARGS}"
      BUILD_BYPRODUCTS "${ARG_BUILD_BYPRODUCTS}"
      CONFIGURE_HANDLED_BY_BUILD ${ARG_CONFIGURE_HANDLED_BY_BUILD}
      UPDATE_DISCONNECTED TRUE
      BUILD_ALWAYS FALSE
      LOG_DOWNLOAD ON
      LOG_UPDATE ON
      LOG_PATCH ON
      LOG_CONFIGURE ON
      LOG_BUILD ON
      LOG_INSTALL ON
      LOG_TEST ON
      LOG_MERGED_STDOUTERR ON
      LOG_OUTPUT_ON_FAILURE ON)
  endif()
endfunction()
# ---------------------------------------------------------------------------------
# cmake-format: on
