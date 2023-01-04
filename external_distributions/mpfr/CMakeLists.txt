# ==============================================================================
# Install the MPFR library required by CGAL for multiprecision numbers.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered mpfr subdirectory")

# Look for an already-installed version of MPFR. If OOMPH_USE_MPFR_FROM is
# defined then we'll look under there. (We have our own FindMPFR.cmake that will
# inspect the OOMPH_USE_MPRF_FROM variable)
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
