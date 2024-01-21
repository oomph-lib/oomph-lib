# cmake-format: off
# ~~~
# Tries to find an installion of GKlib
#
# Once done this will define:
#
# GKLIB_FOUND - BOOL: System has the GKlib library installed
# GKLIB_INCLUDE_DIRS - LIST: The GKlib include directories
# GKLIB_C_LIBRARIES - LIST: The libraries needed to use GKlib via it's C interface
# GKLIB_CXX_LIBRARIES - LIST: The libraries needed to use GKlib via it's C++ interface
# ~~~

include(FindPackageHandleStandardArgs)

# See if the user specified a custom GKLIB installation location using OOMPH_USE_GKLIB_FROM or GKLIBDIR
find_library(GKLIB_C_LIBRARIES NAMES GKlib PATHS "${OOMPH_USE_GKLIB_FROM}/lib" NO_DEFAULT_PATH DOC "GKlib C libraries")
find_path(GKLIB_C_INCLUDES NAMES GKlib.h PATHS "${OOMPH_USE_GKLIB_FROM}/include" NO_DEFAULT_PATH DOC "GKlib C header")

# Try to find libraries and headers in the standard system paths
find_library(GKLIB_C_LIBRARIES NAMES GKlib PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "GKLIB C libraries")
find_path(GKLIB_C_INCLUDES NAMES gmp.h PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "GKLIB C header")

# Handle QUIET and REQUIRED and check the necessary variables were set and if so
# set ``GKLIB_FOUND``
find_package_handle_standard_args(GKlib REQUIRED_VARS GKLIB_C_LIBRARIES
                                                      GKLIB_C_INCLUDES)

if(GKLIB_FOUND)
  set(GKLIB_INCLUDE_DIRS "${GKLIB_C_INCLUDES}")
  list(REMOVE_DUPLICATES GKLIB_INCLUDE_DIRS)

  # Make sure all of the expected headers are in the include directory
  include(CheckIncludeFile)
  foreach(GKLIB_HEADER_FILE IN ITEMS GKlib.h gk_arch.h gk_defs.h gk_externs.h gk_getopt.h gk_macros.h gk_mkblas.h gk_mkmemory.h gk_mkpqueue.h gk_mkpqueue2.h gk_mkrandom.h gk_mksort.h gk_mkutils.h gk_ms_inttypes.h gk_ms_stat.h gk_ms_stdint.h gk_proto.h gk_struct.h gk_types.h gkregex.h)
    find_file(FOUND_HEADER_FILE NAMES ${GKLIB_HEADER_FILE} PATHS ${GKLIB_INCLUDE_DIRS} REQUIRED)
    if(NOT FOUND_HEADER_FILE)
      message(FATAL_ERROR "Header file '${GKLIB_INCLUDE_DIRS}/${GKLIB_HEADER_FILE}' does not exist!")
    endif()
  endforeach()

  if(NOT TARGET GKlib::GKlib)
    add_library(GKlib::GKlib UNKNOWN IMPORTED)
    set_target_properties(
      GKlib::GKlib
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GKLIB_INCLUDE_DIRS}"
                 IMPORTED_LOCATION "${GKLIB_C_LIBRARIES}")
  endif()
endif()
# cmake-format: on
