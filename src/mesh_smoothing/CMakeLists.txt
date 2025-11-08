# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mesh_smoothing subdirectory")

# Define the headers
set(HEADERS mesh_smooth.h)

# Define the sources
set(SOURCES)

# Set the name of the library
set(LIBNAME mesh_smoothing)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving mesh_smoothing subdirectory")
# ------------------------------------------------------------------------------
