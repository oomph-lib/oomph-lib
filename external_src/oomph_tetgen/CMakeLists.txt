# ==============================================================================
# Tetgen installed as library
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oomph_tetgen subdirectory")

# ------------------------------------------------------------------------------
# Retained from Automake build:
#
# There is some hackery going on here because on my main test machine I kept
# getting seg faults when tetgen.cxx was compiled with full optimisation and
# predicates.cxx is supposed to be compiled with no optimisation. Hopefully, one
# day automake will sort this problem out with program specific flags, but it
# may take some time.
#
# The problem may be to do with some strange initialisations within loops and
# loop unrolling, but that could be a complete red herring.
# ------------------------------------------------------------------------------
# Snippet from predicates.cxx:
#
# "...Also be sure to turn on the optimizer when compiling this file."
#
# Hmmm...
# ------------------------------------------------------------------------------

# Define the sources
set(SOURCES tetgen.cxx predicates.cxx)
set_source_files_properties(predicates.cxx PROPERTIES COMPILE_FLAGS -O0)

# Define the headers
set(HEADERS tetgen.h)

# Define the name of the library
set(LIBNAME oomph_tetgen)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

# Silence the annoying warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(
  TARGET ${LIBNAME} CXX_COMPILE_FLAGS -Wno-int-to-pointer-cast
                                      -Wno-unused-result)

# ------------------------------------------------------------------------------
# BOILERPLATE: Update variables in parent scope

# Update list of external sources preprocessor definitions and list of libraries
list(APPEND EXTERNAL_SRC_CXX_DEFINITIONS OOMPH_HAS_TETGEN_LIB)
list(APPEND EXTERNAL_SRC_LIBS ${LIBNAME})

# Make the updated variable values visible to the parent scope
set(EXTERNAL_SRC_CXX_DEFINITIONS ${EXTERNAL_SRC_CXX_DEFINITIONS} PARENT_SCOPE)
set(EXTERNAL_SRC_LIBS ${EXTERNAL_SRC_LIBS} PARENT_SCOPE)

# Define an internal project variable to indicate that we possess this library
set(OOMPH_HAS_TETGEN TRUE CACHE INTERNAL "")

message(VERBOSE "Leaving oomph_tetgen subdirectory")
# ------------------------------------------------------------------------------
