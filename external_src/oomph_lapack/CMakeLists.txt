# ==============================================================================
# Lapack sources required by SuperLU. We turn it into the library oomph_lapack
# to avoid clashes with any existing installations of SuperLU (likely to be
# accessible via -llapack)
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_lapack subdirectory")

# ------------------------------------------------------------------------------
# TODO: Find out whether we need to build statically. From Makefile.am (but
# commented out):
#
# liboomph_lapack_la_LDFLAGS = -static
# ------------------------------------------------------------------------------

# Define the headers
set(HEADERS oomph_required_lapack_routines.h)

# Define the sources
set(SOURCES lsame.c xerbla.c)

# Define the name of the library
set(LIBNAME oomph_lapack)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# Silence the annoying warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(
  TARGET ${LIBNAME} C_COMPILE_FLAGS -Wno-implicit-function-declaration)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Add the library to the list of external source libraries and make the updated
# variable value visible to the parent scope
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

message(VERBOSE "Leaving oomph_lapack subdirectory")
# ------------------------------------------------------------------------------
