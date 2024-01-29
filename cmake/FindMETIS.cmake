# ~~~
# Tries to find an installion of METIS
#
# Once done this will define:
#
# METIS_FOUND - BOOL: System has the METIS library installed
# METIS_INCLUDE_DIRS - LIST: The METIS include directories
# METIS_LIBRARIES - LIST: The libraries needed to use METIS
# ~~~

include(FindPackageHandleStandardArgs)

# cmake-format: off
# See if the user specified a custom METIS installation location using OOMPH_USE_METIS_FROM or METISDIR
find_library(METIS_LIBRARIES NAMES metis PATHS "${OOMPH_USE_METIS_FROM}/lib" NO_DEFAULT_PATH DOC "METIS library")
find_path(METIS_INCLUDES NAMES metis.h PATHS "${OOMPH_USE_METIS_FROM}/include" NO_DEFAULT_PATH DOC "METIS header")

# Try to find libraries and headers in the standard system paths
find_library(METIS_LIBRARIES NAMES metis PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "METIS library")
find_path(METIS_INCLUDES NAMES metis.h PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "METIS header")
# cmake-format: on

# Handle QUIET and REQUIRED and check the necessary variables were set and if so
# set ``METIS_FOUND``
find_package_handle_standard_args(METIS REQUIRED_VARS METIS_LIBRARIES
                                                      METIS_INCLUDES)

if(METIS_FOUND)
  set(METIS_INCLUDE_DIRS "${METIS_INCLUDES}")
  list(REMOVE_DUPLICATES METIS_INCLUDE_DIRS)

  if(NOT TARGET METIS::METIS)
    add_library(METIS::METIS UNKNOWN IMPORTED)
    set_target_properties(
      METIS::METIS
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIRS}"
                 IMPORTED_LOCATION "${METIS_LIBRARIES}")
  endif()
endif()
