# Tries to find an install of the GNU multiple precision library

include(FindPackageHandleStandardArgs)

# cmake-format: off
# See if the user specified a custom GMP installation location using OOMPH_USE_GMP_FROM or GMPDIR
find_library(GMP_C_LIBRARIES NAMES gmp PATHS "${OOMPH_USE_GMP_FROM}/lib" "$ENV{GMPDIR}/lib" NO_DEFAULT_PATH DOC "GMP C libraries")
find_path(GMP_C_INCLUDES NAMES gmp.h PATHS "${OOMPH_USE_GMP_FROM}/include" "$ENV{GMPDIR}/include" NO_DEFAULT_PATH DOC "GMP C header")

# Try to find libraries and headers in the standard system paths
find_library(GMP_C_LIBRARIES NAMES gmp PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "GMP C libraries")
find_path(GMP_C_INCLUDES NAMES gmp.h PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar DOC "GMP C header")
# cmake-format: on

# Handle QUIET and REQUIRED and check the necessary variables were set and if so
# set ``GMP_FOUND``
find_package_handle_standard_args(GMP REQUIRED_VARS GMP_C_LIBRARIES
                                                    GMP_C_INCLUDES)

if(GMP_FOUND)
  set(GMP_INCLUDE_DIRS "${GMP_C_INCLUDES}")
  list(REMOVE_DUPLICATES GMP_INCLUDE_DIRS)

  if(NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    set_target_properties(
      GMP::GMP PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIRS}"
                          IMPORTED_LOCATION "${GMP_C_LIBRARIES}")
  endif()
endif()
