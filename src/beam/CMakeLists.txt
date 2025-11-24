# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered beam subdirectory")

# Define the headers
set(HEADERS beam_elements.h)

# Define the sources
set(SOURCES beam_elements.cc)

# Set the name of the library
set(LIBNAME beam)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving beam subdirectory")
# ------------------------------------------------------------------------------
