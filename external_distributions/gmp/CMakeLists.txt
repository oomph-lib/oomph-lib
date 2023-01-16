# ==============================================================================
# Install the GMP library required by CGAL for multiprecision numbers.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered gmp subdirectory")

# Look for an already-installed version of GMP using FindGMP.cmake. If
# OOMPH_USE_GMP_FROM is defined then FindGMP.cmake will look under there for the
# required files.
find_package(GMP 4.2 REQUIRED)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_GMP)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_GMP TRUE CACHE INTERNAL "")

message(STATUS "Leaving gmp subdirectory")
# ------------------------------------------------------------------------------
