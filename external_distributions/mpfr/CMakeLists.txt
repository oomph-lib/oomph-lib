# ==============================================================================
# Install the MPFR library required by CGAL for multiprecision numbers.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered mpfr subdirectory")

# Look for an already-installed version of MPFR using FindMPFR.cmake. If
# OOMPH_USE_MPFR_FROM is defined then FindMPFR.cmake will look under there for
# the required files.
find_package(MPFR 2.2.1 REQUIRED)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_MPFR)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_MPFR TRUE CACHE INTERNAL "")

message(STATUS "Leaving mpfr subdirectory")
# ------------------------------------------------------------------------------
