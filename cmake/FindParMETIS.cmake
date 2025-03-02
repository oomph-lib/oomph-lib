# ~~~
# Tries to find an installation of ParMETIS
#
# Defines:
#   PARMETIS_FOUND - BOOL: System has the ParMETIS library installed
#   PARMETIS_INCLUDE_DIRS - LIST: The ParMETIS include directories
#   PARMETIS_LIBRARIES - LIST: The libraries needed to use ParMETIS
# ~~~

include(FindPackageHandleStandardArgs)

# Look for user-specified ParMETIS installation
find_library(
  PARMETIS_LIBRARIES
  NAMES parmetis
  PATHS "${OOMPH_USE_PARMETIS_FROM}/lib"
  NO_DEFAULT_PATH
  DOC "PARMETIS libraries")

find_path(
  PARMETIS_INCLUDES
  NAMES parmetis.h
  PATHS "${OOMPH_USE_PARMETIS_FROM}/include"
  NO_DEFAULT_PATH
  DOC "PARMETIS header")

# Search in system paths if not found in custom paths
find_library(
  PARMETIS_LIBRARIES
  NAMES parmetis
  PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar
  DOC "PARMETIS libraries")

find_path(
  PARMETIS_INCLUDES
  NAMES parmetis.h
  PATHS /usr/local /usr /opt/homebrew/opt /usr/local/Cellar
  DOC "PARMETIS header")

# Handle standard package checks
find_package_handle_standard_args(ParMETIS REQUIRED_VARS PARMETIS_LIBRARIES
                                                         PARMETIS_INCLUDES)

if(PARMETIS_FOUND)
  set(PARMETIS_INCLUDE_DIRS "${PARMETIS_INCLUDES}")
  list(REMOVE_DUPLICATES PARMETIS_INCLUDE_DIRS)

  if(NOT TARGET ParMETIS::ParMETIS)
    # Extract the file extension
    get_filename_component(LIBRARY_EXTENSION ${PARMETIS_LIBRARIES} EXT)

    if(LIBRARY_EXTENSION STREQUAL ".a")
      set(LIBRARY_TYPE STATIC)
    elseif(LIBRARY_EXTENSION MATCHES ".(so|dylib)")
      set(LIBRARY_TYPE SHARED)
    else()
      message(
        FATAL_ERROR
          "Unknown library extension for ${PARMETIS_LIBRARIES}. Supported extensions are .a, .so, .dylib"
      )
    endif()

    # Ensure that the SuperLU_DIST target is added only once
    add_library(ParMETIS::ParMETIS ${LIBRARY_TYPE} IMPORTED)
    set_target_properties(
      ParMETIS::ParMETIS
      PROPERTIES IMPORTED_LOCATION "${PARMETIS_LIBRARIES}"
                 INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDES}")
  endif()
endif()
