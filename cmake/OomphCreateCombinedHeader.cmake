# cmake-format: off
# =============================================================================
# Create a combined header from the input headers.
#
# The combined header will contain "#include<...>" commands for the real headers
# in the specified subdirectory. This module is used to emulate the behaviour of
# the include_HEADERS target with Autotools.
#
# Usage:
#
#  include(OomphCreateCombinedHeader)
#  oomph_create_combined_header(TARGET <full-path-of-combined-header>
#                         HEADERS <list-of-headers>
#                         SUBDIRECTORY <subdirectory-of-include-directory>)
# =============================================================================
# cmake-format: on
include_guard()

function(oomph_create_combined_header)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS "")
  set(SINGLE_VALUE_ARGS TARGET SUBDIRECTORY)
  set(MULTI_VALUE_ARGS HEADERS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the prefixed arguments in this scope with clearer variable names
  set(HEADERS ${${PREFIX}_HEADERS})
  set(TARGET ${${PREFIX}_TARGET})
  set(SUBDIRECTORY ${${PREFIX}_SUBDIRECTORY})

  # Construct the header for the combined header file
  set(COMBINED_HEADER_FILE
      "// This file was generated automatically during the CMake process and\n")
  string(APPEND COMBINED_HEADER_FILE "// it will be remade automatically\n")

  # Add the "#include<...>" commands
  foreach(FILE IN LISTS HEADERS)
    # The header path might be a long path (e.g. corresponding to a path in the
    # build tree) but once installed, the header will just be in SUBDIRECTORY,
    # so we have to remember to strip the path before constructing the path to
    # add to the combined header
    cmake_path(GET FILE FILENAME FILE_NAME)
    string(APPEND COMBINED_HEADER_FILE
           "#include <${SUBDIRECTORY}/${FILE_NAME}>\n")
  endforeach()

  # Write to file
  file(WRITE ${TARGET} ${COMBINED_HEADER_FILE})
endfunction()
