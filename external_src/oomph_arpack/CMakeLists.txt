# ==============================================================================
# DOUBLE PRECISION SOURCES FROM ARPACK96 + PATCHES: We turn these into the
# library oomph_arpack to avoid any clashes with existing installations of
# arpack (probably installed as -larpack)
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_arpack subdirectory")

# Define the headers
set(HEADERS_NO_COMBINE debug.h stat.h)

# Define the library headers
set(HEADERS)

# Define the sources
set(SOURCES dummy.cc)

# The location of the private external sources
set(PRIVATE_EXTERNAL_SRC "${CMAKE_SOURCE_DIR}/private/external_src/arpack")

# Helper to search for external sources
include(OomphFindExternalSrcFile)

# Look for the headers. If they're not here, look in private external sources.
# If they're not found, create dummy replacements with the same name
foreach(EXTRA_HEADER ${HEADERS_NO_COMBINE})
  oomph_find_external_src_file(
    FILE ${EXTRA_HEADER}
    SEARCH_PATHS "${PRIVATE_EXTERNAL_SRC}"
    DUMMY_FILE ${EXTRA_HEADER}
    TOUCH_DUMMY_FILE)
endforeach()

# Look for the all_arpack_sources.f file. If it can't be found in the current
# directory, look in the private external sources
oomph_find_external_src_file(
  FILE all_arpack_sources.f
  TARGET SOURCES
  SEARCH_PATHS "${PRIVATE_EXTERNAL_SRC}"
  SET_FLAG HAVE_ARPACK_SOURCES
  DUMMY_FILE dummy_all_arpack_sources.f)

# Set the library name
set(LIBNAME oomph_arpack)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  HEADERS_NO_COMBINE ${HEADERS_NO_COMBINE})

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_arpack subdirectory")
# -----------------------------------------------------------------------------
