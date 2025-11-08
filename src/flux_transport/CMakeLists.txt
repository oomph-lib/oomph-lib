# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered flux_transport subdirectory")

# Define the headers
set(HEADERS flux_transport_elements.h scalar_advection_elements.h
    euler_elements.h)

# Define the sources
set(SOURCES flux_transport_elements.cc scalar_advection_elements.cc
    euler_elements.cc)

# Set the name of the library
set(LIBNAME flux_transport)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving flux_transport subdirectory")
# ------------------------------------------------------------------------------
