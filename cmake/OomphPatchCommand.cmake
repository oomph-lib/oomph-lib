# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Provides a simple patch command; this command is implemented mainly for use
# in calls to FetchContent and ExternalProject.
#
# USAGE:
# ------
#  oomph_git_patch_command(PATCH_DIRECTORY  <directory-containing-patch-file>
#                          PATCH_FILE       <name-of-patch-file>)
#
# EXAMPLE:
# --------
#  set(PATCH_FILE 0001-Patch-of-third-party-library.patch)
#  oomph_git_patch_command(PATCH_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/patches/
#                          PATCH_FILE ${PATCH_FILE})
#
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(oomph_git_patch_command)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS)
  set(SINGLE_VALUE_ARGS PATCH_DIRECTORY PATCH_FILE)
  set(MULTI_VALUE_ARGS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the variables in this scope without a prefix for clarity
  set(PATCH_DIRECTORY ${${PREFIX}_PATCH_DIRECTORY})
  set(PATCH_FILE ${${PREFIX}_PATCH_FILE})

  # The full path to the patch file
  set(PATH_TO_PATCH_FILE "${PATCH_DIRECTORY}/${PATCH_FILE}")

  # Make sure the arguments are valid
  if(NOT IS_DIRECTORY PATCH_DIRECTORY)
    message(FATAL_ERROR "Patch directory does not exist!")
  elseif(NOT EXISTS PATH_TO_PATCH_FILE)
    message(FATAL_ERROR "Patch file does not exist!")
  endif()

  # FIXME: Need to update the "git checkout <SOURCE_DIR>/CMakeLists.txt" line
  # cmake-format: off
  set(PATCH_COMMAND
      ${CMAKE_COMMAND} -E copy_if_different ${PATH_TO_PATCH_FILE} <SOURCE_DIR>
      && git checkout <SOURCE_DIR>/CMakeLists.txt
      && git apply <SOURCE_DIR>/${PATCH_FILE})
  # cmake-format: on
endfunction()
# ------------------------------------------------------------------------------
