# cmake-format: off
# ~~~
# Tries to find an installion of SuperLU_Dist
#
# Once done this will define:
#
# SUPERLU_DIST_FOUND - BOOL: System has the SuperLU_Dist library installed
# SUPERLU_DIST_INCLUDE_DIRS - LIST: The SuperLU_Dist include directories
# SUPERLU_DIST_LIBRARIES - LIST: The libraries needed to use SuperLU_Dist
# ~~~

include(FindPackageHandleStandardArgs)

# See if the user specified a custom SUPERLU_DIST installation location using OOMPH_USE_SUPERLU_DIST_FROM or SUPERLU_DISTDIR
find_library(SUPERLU_DIST_LIBRARIES NAMES superlu_dist PATHS "${OOMPH_USE_SUPERLU_DIST_FROM}/lib" NO_DEFAULT_PATH DOC "SuperLU_DIST C libraries")
find_path(SUPERLU_DIST_INCLUDES NAMES superlu_dist_config.h PATHS "${OOMPH_USE_SUPERLU_DIST_FROM}/include" NO_DEFAULT_PATH DOC "SuperLU_DIST C header")

# Try to find libraries and headers in the standard system paths
find_library(SUPERLU_DIST_LIBRARIES NAMES superlu_dist PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "SuperLU_DIST C libraries")
find_path(SUPERLU_DIST_INCLUDES NAMES superlu_dist_config.h PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "SuperLU_DIST C header")

# Handle QUIET and REQUIRED and check the necessary variables were set and if so
# set ``SUPERLU_DIST_FOUND``
find_package_handle_standard_args(SuperLU_DIST REQUIRED_VARS SUPERLU_DIST_LIBRARIES
                                                      SUPERLU_DIST_INCLUDES)

if(SUPERLU_DIST_FOUND)
  set(SUPERLU_DIST_INCLUDE_DIRS "${SUPERLU_DIST_INCLUDES}")
  list(REMOVE_DUPLICATES SUPERLU_DIST_INCLUDE_DIRS)

  # Make sure all of the expected headers are in the include directory
  include(CheckIncludeFile)
  foreach(SUPERLU_DIST_HEADER_FILE IN ITEMS dcomplex.h machines.h psymbfact.h superlu_FCnames.h superlu_FortranCInterface.h superlu_ddefs.h superlu_defs.h superlu_dist_config.h superlu_enum_consts.h superlu_sdefs.h superlu_zdefs.h supermatrix.h util_dist.h)
    find_file(FOUND_HEADER_FILE NAMES ${SUPERLU_DIST_HEADER_FILE} PATHS ${SUPERLU_DIST_INCLUDE_DIRS} REQUIRED)
    if(NOT FOUND_HEADER_FILE)
      message(FATAL_ERROR "Header file '${SUPERLU_DIST_INCLUDE_DIRS}/${SUPERLU_DIST_HEADER_FILE}' does not exist!")
    endif()
  endforeach()

  if(NOT TARGET SuperLU_DIST::SuperLU_DIST)
    add_library(SuperLU_DIST::SuperLU_DIST UNKNOWN IMPORTED)
    set_target_properties(
      SuperLU_DIST::SuperLU_DIST
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_DIST_INCLUDE_DIRS}"
                 IMPORTED_LOCATION "${SUPERLU_DIST_LIBRARIES}")
  endif()
endif()
# cmake-format: on
