# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rigid_body subdirectory")

# Define the headers
set(HEADERS immersed_rigid_body_elements.h)

# Define the sources
set(SOURCES immersed_rigid_body_elements.cc)

# Set the name of the library
set(LIBNAME rigid_body)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving rigid_body subdirectory")
# ------------------------------------------------------------------------------
