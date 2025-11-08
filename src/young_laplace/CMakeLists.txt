# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered young_laplace subdirectory")

# Define the headers
set(HEADERS young_laplace_elements.h young_laplace_contact_angle_elements.h
    refineable_young_laplace_elements.h)

# Define the sources
set(SOURCES young_laplace_elements.cc young_laplace_contact_angle_elements.cc
    refineable_young_laplace_elements.cc)

# Set the name of the library
set(LIBNAME young_laplace)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving young_laplace subdirectory")
# ------------------------------------------------------------------------------
