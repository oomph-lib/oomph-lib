# --------------------------------------------------------------------------------
# Try to find MPFR
#
# Once done this will define:
#
# * MPFR_FOUND        - System has MPFR
# * MPFR_INCLUDES_DIRS - The MPFR include directories
# * MPFR_LIBRARIES    - The libraries needed to use MPFR
# --------------------------------------------------------------------------------

# cmake-format: off
# See if the user specified a custom MPFR installation location using OOMPH_USE_MPFR_FROM
# or MPFRDIR. If not try to find the library/headers in the standard system paths
find_library(MPFR_LIBRARIES mpfr PATHS "${OOMPH_USE_MPFR_FROM}/lib" "$ENV{MPFRDIR}/lib" NO_DEFAULT_PATH DOC "MPFR C library")
find_library(MPFR_LIBRARIES mpfr PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "MPFR C library")
find_path(MPFR_INCLUDES NAMES mpfr.h PATHS "${OOMPH_USE_MPFR_FROM}/include" "$ENV{MPFRDIR}/include" NO_DEFAULT_PATH DOC "MPFR C header")
find_path(MPFR_INCLUDES NAMES mpfr.h PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "MPFR C header")
# cmake-format: on

# Handle the QUIET and REQUIRED arguments and set MPFR_FOUND to TRUE if all
# listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MPFR
  FOUND_VAR MPFR_FOUND
  REQUIRED_VARS MPFR_LIBRARIES MPFR_INCLUDES
  FAIL_MESSAGE "Could NOT find MPFR, use MPFR_ROOT to hint its location")

if(MPFR_FOUND)
  set(MPFR_INCLUDE_DIRS "${MPFR_INCLUDES}")
  list(REMOVE_DUPLICATES MPFR_INCLUDE_DIRS)

  if(NOT TARGET MPFR::MPFR)
    add_library(MPFR::MPFR UNKNOWN IMPORTED)
    set_target_properties(
      MPFR::MPFR PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDES}"
                            IMPORTED_LOCATION "${MPFR_LIBRARIES}")
  endif()
endif()
