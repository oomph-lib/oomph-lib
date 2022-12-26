# ==============================================================================
# Install the MPFR library required by CGAL for multiprecision numbers.
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered mpfr subdirectory")

# Look for an already-installed version of MPFR. If OOMPH_USE_MPFR_FROM is
# defined then we'll look under there
if(OOMPH_USE_MPFR_FROM)
  find_package(
    MPFR 3.1.6 REQUIRED
    PATHS "${OOMPH_USE_MPFR_FROM}"
    NO_DEFAULT_PATH)
endif()
find_package(MPFR 3.1.6 REQUIRED)

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
