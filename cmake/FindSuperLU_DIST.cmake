# cmake-format: off
# ------------------------------------------------------------------------------
# Tries to find an installation of SuperLU_Dist
#
# Once done, this will define:
#   SUPERLU_DIST_FOUND - BOOL: System has the SuperLU_Dist library installed
#   SUPERLU_DIST_INCLUDE_DIRS - LIST: The SuperLU_Dist include directories
#   SUPERLU_DIST_LIBRARIES - LIST: The libraries needed to use SuperLU_Dist
#
# ------------------------------------------------------------------------------
# cmake-format: on
include(FindPackageHandleStandardArgs)

# Define search paths for libraries and headers, adding common locations for
# both macOS and Linux
find_library(
  SUPERLU_DIST_LIBRARIES
  NAMES superlu_dist
  PATHS "${OOMPH_USE_SUPERLU_DIST_FROM}/lib" /usr/local /usr /opt/homebrew/opt
        /usr/local/Cellar
  NO_DEFAULT_PATH
  DOC "SuperLU_DIST C libraries")

find_path(
  SUPERLU_DIST_INCLUDES
  NAMES superlu_dist_config.h
  PATHS "${OOMPH_USE_SUPERLU_DIST_FROM}/include" /usr/local /usr
        /opt/homebrew/opt /usr/local/Cellar
  DOC "SuperLU_DIST C header")

# Handle QUIET and REQUIRED and check the necessary variables were set and if so
# set `SUPERLU_DIST_FOUND`
find_package_handle_standard_args(
  SuperLU_DIST REQUIRED_VARS SUPERLU_DIST_LIBRARIES SUPERLU_DIST_INCLUDES)

if(SUPERLU_DIST_FOUND)
  set(SUPERLU_DIST_INCLUDE_DIRS "${SUPERLU_DIST_INCLUDES}")
  list(REMOVE_DUPLICATES SUPERLU_DIST_INCLUDE_DIRS)

  # Make sure all of the expected headers are in the include directory
  include(CheckIncludeFile)
  foreach(SUPERLU_DIST_HEADER_FILE IN ITEMS superlu_ddefs.h
                                            superlu_enum_consts.h)
    find_file(
      FOUND_HEADER_FILE
      NAMES ${SUPERLU_DIST_HEADER_FILE}
      PATHS ${SUPERLU_DIST_INCLUDE_DIRS} REQUIRED)
    if(NOT FOUND_HEADER_FILE)
      message(
        FATAL_ERROR
          "Header file '${SUPERLU_DIST_INCLUDE_DIRS}/${SUPERLU_DIST_HEADER_FILE}' does not exist!"
      )
    endif()
  endforeach()

  # Ensure that the SuperLU_DIST target is added only once
  if(NOT TARGET SuperLU_DIST::SuperLU_DIST)
    # Extract the file extension from the SUPERLU_DIST_LIBRARIES path
    get_filename_component(LIBRARY_EXTENSION ${SUPERLU_DIST_LIBRARIES} EXT)

    if(LIBRARY_EXTENSION STREQUAL ".a")
      set(LIBRARY_TYPE STATIC)
    elseif(LIBRARY_EXTENSION MATCHES ".(so|dylib)")
      set(LIBRARY_TYPE SHARED)
    else()
      message(
        FATAL_ERROR
          "Unknown library extension for ${SUPERLU_DIST_LIBRARIES}. Supported extensions are .a, .so, .dylib"
      )
    endif()

    # Ensure that the SuperLU_DIST target is added only once
    add_library(SuperLU_DIST::SuperLU_DIST ${LIBRARY_TYPE} IMPORTED)
    set_target_properties(
      SuperLU_DIST::SuperLU_DIST
      PROPERTIES IMPORTED_LOCATION "${SUPERLU_DIST_LIBRARIES}"
                 INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_DIST_INCLUDE_DIRS}")
  endif()
endif()
