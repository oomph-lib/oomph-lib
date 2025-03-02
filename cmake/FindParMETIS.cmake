# ~~~
# Tries to find an installion of ParMETIS
#
# Once done this will define:
#
# PARMETIS_FOUND - BOOL: System has the ParMETIS library installed
# PARMETIS_INCLUDE_DIRS - LIST: The ParMETIS include directories
# PARMETIS_LIBRARIES - LIST: The libraries needed to use ParMETIS
# ~~~

include(FindPackageHandleStandardArgs)

# cmake-format: off
# See if the user specified a custom PARMETIS installation location using OOMPH_USE_PARMETIS_FROM or PARMETISDIR
find_library(PARMETIS_LIBRARIES NAMES parmetis PATHS "${OOMPH_USE_PARMETIS_FROM}/lib" NO_DEFAULT_PATH DOC "PARMETIS libraries")
find_path(PARMETIS_INCLUDES NAMES parmetis.h PATHS "${OOMPH_USE_PARMETIS_FROM}/include" NO_DEFAULT_PATH DOC "PARMETIS header")

# Try to find libraries and headers in the standard system paths
find_library(PARMETIS_LIBRARIES NAMES parmetis PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "PARMETIS libraries")
find_path(PARMETIS_INCLUDES NAMES parmetis.h PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "PARMETIS header")
# cmake-format: on

# Handle QUIET and REQUIRED and check the necessary variables were set and if so
# set ``PARMETIS_FOUND``
find_package_handle_standard_args(ParMETIS REQUIRED_VARS PARMETIS_LIBRARIES
                                                         PARMETIS_INCLUDES)

if(PARMETIS_FOUND)
  set(PARMETIS_INCLUDE_DIRS "${PARMETIS_INCLUDES}")
  list(REMOVE_DUPLICATES PARMETIS_INCLUDE_DIRS)

  if(NOT TARGET ParMETIS::ParMETIS)
    add_library(ParMETIS::ParMETIS UNKNOWN IMPORTED)
    set_target_properties(
      ParMETIS::ParMETIS
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIRS}"
                 IMPORTED_LOCATION "${PARMETIS_LIBRARIES}")
  endif()
endif()
