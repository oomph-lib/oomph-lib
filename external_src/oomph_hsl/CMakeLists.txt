# ==============================================================================
# Frontal solver from the HSL library We turn it into the library oomph_hsl to
# avoid clashes with other installations of the hsl library, likely to be
# accessible via -lhsl
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_hsl subdirectory")

# Define the headers
set(HEADERS)

# Define the sources
set(SOURCES dummy.cc)

# Look for the frontal.f file. If it can't be found in the current directory,
# look in the private external sources
include(OomphFindExternalSrcFile)
oomph_find_external_src_file(
  FILE frontal.f
  TARGET SOURCES
  SEARCH_PATHS ${CMAKE_SOURCE_DIR}/private/external_src/hsl
  SET_FLAG HAVE_HSL_SOURCES
  DUMMY_FILE dummy_frontal.f)

# Define the name of the library
set(LIBNAME oomph_hsl)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_hsl subdirectory")
# -----------------------------------------------------------------------------
