# cmake-format: off
# ==============================================================================
# Build the MUMPS (and possibly Scalapack) library.
#
# Details on the available CMake build options can be found here:
#                   https://github.com/scivision/mumps
# Note, however, that we do not download and build from the repository linked
# above. It just happens to be where helpful build instructions can be found.
#
# Notable build options:
#
#    -Dautobuild=true (default: true): Will download and build a local copy of
#                                 Lapack and/or Scalapack if missing or broken.
#
#    -Dopenmp=true (default: false): Use OpenMP. Note that OpenMP can make MUMPS
#                                    over 10x slower in certain situations. Try
#                                    with and without OpenMP to see which is
#                                    faster for your situation.
# ==============================================================================
# cmake-format: on
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(STATUS "Entered mumps_and_scalapack subdirectory")

# Set MUMPS build options
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
# set(BUILD_TESTING OFF)

# TODO: Find out how to control the installation destination of the stupid thing

# Download MUMPS repo. Currently using a personal forked repo with patched
# changes to address the overriden install path when the default one is used. If
# it ever gets fixed in the future, we should fork a copy of the repo here:
#
# https://github.com/scivision/mumps.git
#
# to the oomph-lib organisation.
include(FetchContent)
FetchContent_Declare(
  mumps_project
  GIT_REPOSITORY https://github.com/scivision/mumps.git
  GIT_TAG v5.4.1.2
  # GIT_REPOSITORY https://github.com/PuneetMatharu/mumps.git GIT_TAG v5.4.0.6
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(mumps_project)

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
