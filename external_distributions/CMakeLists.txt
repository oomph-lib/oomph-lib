# cmake-format: off
# ==============================================================================
#                                          888             888 d8b 888
#                                          888             888 Y8P 888
#                                          888             888     888
#  .d88b.   .d88b.  88888b.d88b.  88888b.  88888b.         888 888 88888b.
# d88""88b d88""88b 888 "888 "88b 888 "88b 888 "88b        888 888 888 "88b
# 888  888 888  888 888  888  888 888  888 888  888 888888 888 888 888  888
# Y88..88P Y88..88P 888  888  888 888 d88P 888  888        888 888 888 d88P
#  "Y88P"   "Y88P"  888  888  888 88888P"  888  888        888 888 88888P"
#                                 888
#                                 888
#                                 888
# ------------------------------------------------------------------------------
# oomph-lib is an object-oriented, open-source finite-element library for the
# simulation of multi-physics problems. It is developed and maintained by
# Matthias Heil and Andrew Hazel of the Department of Mathematics at The
# University of Manchester, along with many other contributors.
# ==============================================================================
# cmake-format: on

# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(
  oomph_lib_third_party_libraries
  DESCRIPTION "Builds the third-party libraries used by oomph-lib."
  VERSION 1.0.0
  HOMEPAGE_URL "https://oomph-lib.github.io/oomph-lib/doc/html/"
  LANGUAGES C CXX Fortran)

# ---------------------------[ INITIAL PROJECT SETUP ]--------------------------

# There is no 'install' step here; just configure and build!
set(CMAKE_SKIP_INSTALL_RULES ON)

# Handle warning about timestamps when downloading archives via ExternalProject
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# Specify the location of non-standard CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/")

# Specify the C++ standard. NOTE: This version MUST align with the version used
# in the build of the main oomph-lib library
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Build in Release (fully optimised) mode by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build mode for performance."
      FORCE)
endif()

# Select level of verbosity (e.g. STATUS, VERBOSE, etc.)
set(CMAKE_MESSAGE_LOG_LEVEL VERBOSE)

# ------------------------------[ BUILD OPTIONS ]-------------------------------
# cmake-format: off
include(OomphCustomOptions)

# Boolean options
oomph_option(OOMPH_ENABLE_MPI "Build third-party libraries with MPI support?" OFF)
oomph_option(OOMPH_BUILD_OPENBLAS "Build OpenBLAS?" ON)
oomph_option(OOMPH_BUILD_SUPERLU "Build SuperLU?" ON)
oomph_option(OOMPH_BUILD_SUPERLU_DIST "Build SuperLUDist?" ${OOMPH_ENABLE_MPI})
oomph_option(OOMPH_BUILD_CGAL "Build CGAL?" ON)
oomph_option(OOMPH_BUILD_MUMPS "Build MUMPS?" ON)
oomph_option(OOMPH_BUILD_HYPRE "Build Hypre?" ON)
oomph_option(OOMPH_BUILD_TRILINOS "Build Trilinos?" ON)
oomph_option(OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS "Disable testing when building the third-party libraries?" OFF)

# Set the default installation path
if(PROJECT_IS_TOP_LEVEL AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/install" CACHE PATH "Install path prefix" FORCE)
endif()

# Path-like flags
oomph_path_option(
  FLAG OOMPH_USE_OPENBLAS_FROM
  DOCSTRING "The path to a preinstalled version of OpenBLAS."
)
oomph_path_option(
  FLAG OOMPH_USE_GKLIB_FROM
  DOCSTRING "The path to a preinstalled version of GKLIB."
)
oomph_path_option(
  FLAG OOMPH_USE_METIS_FROM
  DOCSTRING "The path to a preinstalled version of METIS."
)
oomph_path_option(
  FLAG OOMPH_USE_PARMETIS_FROM
  DOCSTRING "The path to a preinstalled version of ParMETIS."
)
oomph_path_option(
  FLAG OOMPH_USE_BOOST_FROM
  DOCSTRING "The path to a preinstalled version of Boost."
)

# cmake-format: on
# ----------------------[ FIND REQUIRED/DESIRED PACKAGES ]----------------------

# Use CMake built-in modules to search for certain packages
find_package(Git)
find_package(Python3 REQUIRED)

# Helpful CMake modules
include(ExternalProject)
include(ProcessorCount)

# The 'make' executable
find_program(MAKE_EXECUTABLE NAMES make REQUIRED)

# Get our helper function for getting external libraries
include(OomphGetExternalProjectHelper)

# --------------------------------[ GENERAL ]----------------------------------

# Decide on number of jobs to use for building projects
ProcessorCount(OOMPH_NUM_JOBS)
if(OOMPH_NUM_JOBS EQUAL 0)
  set(OOMPH_NUM_JOBS 1)
  message(WARNING "Could not determine CPU cores. Using single job.")
else()
  message(STATUS "Number of jobs for parallel make: ${OOMPH_NUM_JOBS}")
endif()

# ----------------------------------[ MPI ]------------------------------------

# Set up MPI functionality (see cmake/OomphMPI.cmake)
if(OOMPH_ENABLE_MPI)
  include(OomphMPI)
endif()

# --------------------------------[ OPENBLAS ]---------------------------------

if(OOMPH_USE_OPENBLAS_FROM)
  find_package(OpenBLAS 0.3.25 REQUIRED PATHS ${OOMPH_USE_OPENBLAS_FROM})
  set(OOMPH_BUILD_OPENBLAS OFF CACHE PATH "Build OpenBLAS?" FORCE)
elseif(OOMPH_BUILD_OPENBLAS)
  include(OomphGetExternalOpenBLAS)
endif()

# Get the OpenBLAS root directory (required by MUMPS)
if(OpenBLAS_LIBRARIES)
  cmake_path(GET OpenBLAS_LIBRARIES PARENT_PATH OpenBLAS_LIBDIR)
  cmake_path(GET OpenBLAS_LIBDIR PARENT_PATH OpenBLAS_ROOT)
  set(OpenBLAS_ROOT "${OpenBLAS_ROOT}" CACHE PATH "" FORCE)
endif()

if(NOT OpenBLAS_LIBRARIES)
  if(OOMPH_BUILD_HYPRE OR OOMPH_BUILD_TRILINOS)
    message(FATAL_ERROR "Refusing to build HYPRE or Trilinos without OpenBLAS!")
  endif()
endif()

# -------------------------------[ SUPERLUDIST ]-------------------------------

# TODO: Try updating SuperLU to v9

# If we're given SuperLUDist or need to build it
if(OOMPH_BUILD_SUPERLU OR OOMPH_BUILD_SUPERLU_DIST)
  if((OOMPH_BUILD_SUPERLU_DIST) AND (NOT OOMPH_ENABLE_MPI))
    set(OOMPH_BUILD_SUPERLU_DIST "OFF")
    message(
      WARNING
        "Not building SuperLUDist (despite request) because MPI is disabled!")
  endif()

  # If we're given GKlib
  if(OOMPH_USE_GKLIB_FROM)
    find_package(GKlib REQUIRED)
  endif()

  # If we're given METIS
  if(OOMPH_USE_METIS_FROM)
    find_package(METIS REQUIRED)
  endif()

  # If we're given ParMETIS
  if(OOMPH_USE_PARMETIS_FROM)
    find_package(ParMETIS REQUIRED)
  endif()

  # Build SuperLU/SuperLUDist
  include(OomphGetExternalSuperLUAndOrSuperLUDist)
  set(GKLIB_FOUND TRUE)
  set(METIS_FOUND TRUE)
  set(PARMETIS_FOUND TRUE)
  set(SUPERLU_FOUND TRUE)
  set(SUPERLU_DIST_FOUND TRUE)
else()
  # If we haven't been asked to build SuperLUDist but the user (bizarrely) told
  # us where to find GKlib, METIS and/or ParMETIS...
  set(OOMPH_PASSED_SUPERLU_DEPS)
  if(OOMPH_USE_GKLIB_FROM)
    list(APPEND OOMPH_PASSED_SUPERLU_DEPS "'OOMPH_USE_GKLIB_FROM'")
  endif()
  if(OOMPH_USE_METIS_FROM)
    list(APPEND OOMPH_PASSED_SUPERLU_DEPS "'OOMPH_USE_METIS_FROM'")
  endif()
  if(OOMPH_USE_PARMETIS_FROM)
    list(APPEND OOMPH_PASSED_SUPERLU_DEPS "'OOMPH_USE_PARMETIS_FROM'")
  endif()
  if(OOMPH_PASSED_SUPERLU_DEPS)
    string(JOIN " and " OOMPH_PASSED_SUPERLU_DEPS_AS_STR
           ${OOMPH_PASSED_SUPERLU_DEPS})
    message(
      WARNING
        "\nYou specified ${OOMPH_PASSED_SUPERLU_DEPS_AS_STR} but did not enable the build of SuperLU/SuperLUDist. Did you forget to pass '-DOOMPH_BUILD_SUPERLU=ON' or '-DOOMPH_BUILD_SUPERLU_DIST=ON'?"
    )
  endif()
endif()

# ----------------------------------[ CGAL ]-----------------------------------

# If we're given CGAL or need to build it
if(OOMPH_BUILD_CGAL)
  # If we're given Boost. This form of find_package(...) will search for the
  # BoostConfig.cmake shipped with Boost
  if(OOMPH_USE_BOOST_FROM)
    find_package(
      Boost 1.83.0 REQUIRED
      PATHS ${OOMPH_USE_BOOST_FROM}
      NO_DEFAULT_PATH)
  endif()

  # Get CGAL
  include(OomphGetExternalCGAL)
  set(BOOST_FOUND TRUE)
  set(CGAL_FOUND TRUE)
elseif(OOMPH_USE_BOOST_FROM)
  # ...then why did the user tell us where to find Boost?
  message(
    WARNING
      "\nYou specified 'OOMPH_USE_BOOST_FROM' but did not enable the build of CGAL. Did you forget to pass '-DOOMPH_BUILD_CGAL=ON'?"
  )
endif()

# ----------------------------------[ MUMPS ]----------------------------------

if(OOMPH_BUILD_MUMPS)
  include(OomphGetExternalMUMPS)
endif()

# ----------------------------------[ HYPRE ]----------------------------------

if(OOMPH_BUILD_HYPRE)
  include(OomphGetExternalHYPRE)
endif()

# ---------------------------------[ TRILINOS ]--------------------------------

if(OOMPH_BUILD_TRILINOS)
  include(OomphGetExternalTrilinos)
endif()

# ----------------------------[ DOC CONFIGURATION ]----------------------------

include(OomphPrintProjectConfiguration)
oomph_print_project_configuration()

# -------------------------[ HOW TO USE THE LIBRARIES ]------------------------

include(OomphPrintThirdPartyLibrariesUsage)
oomph_print_third_party_libraries_usage()

# -----------------------------------------------------------------------------
