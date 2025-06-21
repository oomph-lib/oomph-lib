# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered womersley subdirectory")

# Define the headers
set(HEADERS womersley_elements.h Twomersley_elements.h)

# Define the sources
set(SOURCES womersley_elements.cc Twomersley_elements.cc)

# Set the name of the library
set(LIBNAME womersley)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving womersley subdirectory")
# ------------------------------------------------------------------------------
