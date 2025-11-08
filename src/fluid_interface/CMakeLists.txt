# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fluid_interface subdirectory")

# Define the headers
set(HEADERS interface_elements.h constrained_volume_elements.h
    specific_node_update_interface_elements.h surfactant_transport_elements.h)

# Define the sources
set(SOURCES interface_elements.cc constrained_volume_elements.cc
    surfactant_transport_elements.cc)

# Set the name of the library
set(LIBNAME fluid_interface)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving fluid_interface subdirectory")
# ------------------------------------------------------------------------------
