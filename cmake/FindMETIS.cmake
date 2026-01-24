# cmake-format: off
# ------------------------------------------------------------------------------
# Tries to find an installation of METIS
#
# Once done, this will define:
#
# METIS_FOUND - BOOL: System has the METIS library installed
# METIS_INCLUDE_DIRS - LIST: The METIS include directories
# METIS_LIBRARIES - LIST: The libraries needed to use METIS
# METIS_VERSION - STRING: Version of METIS (if found)
# ------------------------------------------------------------------------------
# cmake-format: on
include(FindPackageHandleStandardArgs)

# Define common search paths
set(METIS_SEARCH_PATHS
    ${OOMPH_USE_METIS_FROM}
    $ENV{METISDIR}
    /usr/local
    /usr
    /opt/homebrew/opt
    /usr/local/Cellar
    /opt/local)

# Search for the METIS library
find_library(
  METIS_LIBRARIES
  NAMES metis
  PATHS ${METIS_SEARCH_PATHS}
  PATH_SUFFIXES lib lib64
  DOC "METIS library")

# Search for the METIS headers
find_path(
  METIS_INCLUDE_DIRS
  NAMES metis.h
  PATHS ${METIS_SEARCH_PATHS}
  PATH_SUFFIXES include
  DOC "METIS header directory")

# Check for a version file
if(METIS_INCLUDE_DIRS AND EXISTS "${METIS_INCLUDE_DIRS}/metis.h")
  file(STRINGS "${METIS_INCLUDE_DIRS}/metis.h" METIS_VERSION_LINE
       REGEX "#define METIS_VER_MAJOR.*")
  string(REGEX MATCH "[0-9]+" METIS_VERSION_MAJOR "${METIS_VERSION_LINE}")

  file(STRINGS "${METIS_INCLUDE_DIRS}/metis.h" METIS_VERSION_LINE
       REGEX "#define METIS_VER_MINOR.*")
  string(REGEX MATCH "[0-9]+" METIS_VERSION_MINOR "${METIS_VERSION_LINE}")

  set(METIS_VERSION "${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}")
endif()

# Handle REQUIRED and QUIET options
find_package_handle_standard_args(
  METIS
  REQUIRED_VARS METIS_LIBRARIES METIS_INCLUDE_DIRS
  VERSION_VAR METIS_VERSION)

# If METIS was found, create an imported target
if(METIS_FOUND)
  if(NOT TARGET METIS::METIS_HEADER_ONLY)
    # We'll define an INTERFACE version of the library just to propagate the
    # location of metis.h
    add_library(METIS::METIS_HEADER_ONLY INTERFACE IMPORTED)
    set_target_properties(
      METIS::METIS_HEADER_ONLY PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                          "${METIS_INCLUDE_DIRS}")
  endif()

  if(NOT TARGET METIS::METIS)
    # ...and a regular one so that the user can link against the metis lib. It's
    # unlikely we'll need this since the METIS library is only needed by SuperLU
    # and that's built outside of the oomph-lib framework
    add_library(METIS::METIS UNKNOWN IMPORTED)
    set_target_properties(
      METIS::METIS
      PROPERTIES IMPORTED_LOCATION "${METIS_LIBRARIES}"
                 INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIRS}")
  endif()
else()
  message(FATAL_ERROR "Could not find METIS")
endif()
# ------------------------------------------------------------------------------
