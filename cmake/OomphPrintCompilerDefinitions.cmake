# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Prints an alphabetically-sorted list of the compile definitions visible to the
# current scope.
#
# USAGE:
# ------
#          include(OomphPrintCompilerDefinitions)
# =============================================================================
# cmake-format: on

# Get the list compile definitions
get_directory_property(OOMPH_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR}
                                                           COMPILE_DEFINITIONS)

# Sort the list alphabetically
list(SORT OOMPH_COMPILE_DEFINITIONS)

# Now print the compile definitions
message(STATUS "COMPILER DEFINITIONS:")
foreach(DEFN ${OOMPH_COMPILE_DEFINITIONS})
  message(STATUS "  ...${DEFN}")
endforeach()

# ------------------------------------------------------------------------------
