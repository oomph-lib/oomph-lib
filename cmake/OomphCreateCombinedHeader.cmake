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

  # Construct the header for the combined header file
  set(${PREFIX}_COMBINED_HEADER_FILE
      "// This file was generated automatically during the CMake process and\n")
  string(APPEND ${PREFIX}_COMBINED_HEADER_FILE
         "// it will be remade automatically\n")

  # Add the "#include<...>" commands
  foreach(${PREFIX}_FILE ${${PREFIX}_HEADERS})
    string(APPEND ${PREFIX}_COMBINED_HEADER_FILE
           "#include <${${PREFIX}_SUBDIRECTORY}/${${PREFIX}_FILE}>\n")
  endforeach()

  # Write to file
  file(WRITE ${${PREFIX}_TARGET} ${${PREFIX}_COMBINED_HEADER_FILE})
endfunction()
