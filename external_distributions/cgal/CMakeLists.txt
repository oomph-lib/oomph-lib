# ==============================================================================
# Build/install the CGAL library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered cgal subdirectory")

# Make sure the CGAL dependencies have been found
if(NOT OOMPH_HAS_GMP)
  message(FATAL_ERROR "Cannot use CGAL -- the 'gmp' library is missing!")
endif()
if(NOT OOMPH_HAS_MPFR)
  message(FATAL_ERROR "Cannot use CGAL -- the 'mpfr' library is missing!")
endif()
if(NOT OOMPH_HAS_BOOST)
  message(FATAL_ERROR "Cannot use CGAL -- the 'Boost' library is missing!")
endif()

if(OOMPH_USE_CGAL_FROM)
  # If the user specified OOMPH_USE_CGAL_FROM, look there
  find_package(
    CGAL 5.6 REQUIRED
    PATHS "${OOMPH_USE_CGAL_FROM}"
    NO_DEFAULT_PATH)
else()
  # Look in the system paths
  find_package(CGAL 5.6 REQUIRED)
endif()

if(NOT TARGET CGAL::CGAL)
  message(
    FATAL_ERROR
      "Something went wrong during the configuration of the CGAL library!")
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_CGAL)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_CGAL TRUE CACHE INTERNAL "")

message(STATUS "Leaving cgal subdirectory")
# ------------------------------------------------------------------------------
