# cmake-format: off
#=============================================================================
# Try to find the Triangle library:
#
# Once done this will define:
#  TRIANGLE_FOUND   - system has TRIANGLE
#  TRIANGLE_INC     - TRIANGLE include directory
#  TRIANGLE_LIB     - TRIANGLE library
#
# Usage:
#  find_package(TRIANGLE)
#
# Setting these changes the behavior of the search
#  TRIANGLE_FOUND   - Whether the system has Triangle
#  TRIANGLE_INC     - Triangle include directory (location of triangle.h)
#  TRIANGLE_LIB     - Triangle library (libraries to link against to use Triangle)
#=============================================================================
# cmake-format: on

# =============================================================================
# Look for environment variables
# =============================================================================
# Find the location of triangle.h
find_path(TRIANGLE_INC triangle.h DOC "The Triangle include directory")

# Find the Triangle lib files
set(TRIANGLE_NAMES ${TRIANGLE_NAMES} libtriangle triangle)
find_library(
  TRIANGLE_LIB
  NAMES ${TRIANGLE_NAMES}
  DOC "The Triangle library")

# =============================================================================
# CMake check and done
# =============================================================================
# handle the QUIETLY and REQUIRED arguments and set TRIANGLE_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  TRIANGLE
  "Triangle could not be found; make sure that you actually have Triangle!"
  REQUIRED_VARS TRIANGLE_LIB TRIANGLE_INC
  VERSION_VAR TRIANGLE_VERSION)

# if(TRIANGLE_FOUND) set(TRIANGLE_LIB ${TRIANGLE_LIB}) set(TRIANGLE_INC
# ${TRIANGLE_INC}) endif()

message("TRIANGLE LIB/INC: ${TRIANGLE_LIB} ${TRIANGLE_INC}")

mark_as_advanced(TRIANGLE_INC TRIANGLE_LIB)
