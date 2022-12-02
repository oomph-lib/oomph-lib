# ==============================================================================
# Blas sources
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_blas subdirectory")

# Define the headers
set(HEADERS)

# Define the sources
set(SOURCES)

# Define the name of the library
set(LIBNAME oomph_blas)

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

message(VERBOSE "Leaving oomph_blas subdirectory")
# ------------------------------------------------------------------------------
