# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Locate the Gmsh command
#
# USAGE:
# ------
#     include(OomphFindGmsh)
#
# =============================================================================
# cmake-format: on
include_guard()

# ----------------------------------------------------------------------------
if(NOT OOMPH_USE_GMSH_FROM)
  message(
    FATAL_ERROR
      "Tried to search for Gmsh but 'OOMPH_USE_GMSH_FROM' has not been set!")
endif()

# If OOMPH_USE_GMSH_FROM is the path to an executable, make sure it is a
# program. If it's a directory, check that we can find 'gmsh' in the directory
if(NOT EXISTS "${OOMPH_USE_GMSH_FROM}")
  message(
    FATAL_ERROR
      "You said that 'gmsh' can be found here:\n\t${OOMPH_USE_GMSH_FROM}\nbut this file/directory does not exist!"
  )
elseif(IS_DIRECTORY "${OOMPH_USE_GMSH_FROM}")
  message(
    FATAL_ERROR
      "Argument to 'OOMPH_USE_GMSH_FROM' must be an executable, not a directory!"
  )
endif()

# Extract path components
set(OOMPH_GMSH_DIR ${OOMPH_USE_GMSH_FROM})
set(OOMPH_GMSH_PROGRAM)
cmake_path(GET OOMPH_USE_GMSH_FROM STEM OOMPH_GMSH_PROGRAM)
cmake_path(GET OOMPH_USE_GMSH_FROM PARENT_PATH OOMPH_GMSH_DIR)

# Search for a program;
#
# TODO: Add a 'VALIDATOR' argument to check that gmsh works
find_program(
  OOMPH_GMSH_COMMAND
  NAMES gmsh ${OOMPH_GMSH_PROGRAM}
  PATHS ${OOMPH_GMSH_DIR} REQUIRED
  NO_DEFAULT_PATH)

# Notify user
message(STATUS "Command 'gmsh' exists at: ${OOMPH_GMSH_COMMAND}")
# ----------------------------------------------------------------------------
