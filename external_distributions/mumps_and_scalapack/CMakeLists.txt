# cmake-format: off
# ==============================================================================
# Build the MUMPS (and possibly Scalapack) library.
#
# Details on the available CMake build options can be found here:
#               https://github.com/puneetmatharu/mumps.git
# Note, however, that we do not download and build from the repository linked
# above. It just happens to be where helpful build instructions can be found.
# ==============================================================================
# cmake-format: on
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered mumps_and_scalapack subdirectory")

# Set MUMPS configuration options
set(autobuild ON)
set(lapack_external OFF)
set(scalapack_external OFF)
set(mumps_external OFF)
set(parallel ${OOMPH_ENABLE_MPI})
set(intsize64 OFF)
set(metis OFF)
set(scotch OFF)
set(openmp OFF)
set(mumps_matlab OFF)
set(use_default_install_location ON)

# Look for an already-installed version of MUMPS. If OOMPH_USE_MUMPS_FROM is
# defined then we'll look under there
if(OOMPH_USE_MUMPS_FROM)
  find_package(
    MUMPS 5.5.1.11 REQUIRED
    PATHS "${OOMPH_USE_MUMPS_FROM}"
    NO_DEFAULT_PATH)
endif()
find_package(MUMPS 5.5.1.11 QUIET)

# If we couldn't find a pre-installed version, lets build it ourselves
if(MUMPS_FOUND)
  message(STATUS "Found an existing installation of MUMPS at: ${MUMPS_DIR}")
else()
  include(FetchContent)
  FetchContent_Declare(
    mumps_project
    GIT_REPOSITORY https://github.com/puneetmatharu/mumps.git
    # GIT_REPOSITORY https://github.com/scivision/mumps.git
    GIT_TAG v5.5.1.11
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(mumps_project)
endif()

# ------------------------------------------------------------------------------
# BOILERPLATE: Update global variables

# Update external distributions preprocessor definitions and make the updated
# variable values visible to the parent scope
list(APPEND EXTERNAL_DIST_CXX_DEFINITIONS OOMPH_HAS_MUMPS)
set(EXTERNAL_DIST_CXX_DEFINITIONS ${EXTERNAL_DIST_CXX_DEFINITIONS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_MUMPS TRUE CACHE INTERNAL "")

message(STATUS "Leaving mumps_and_scalapack subdirectory")
# ------------------------------------------------------------------------------
