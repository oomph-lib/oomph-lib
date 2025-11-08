# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spherical_navier_stokes subdirectory")

# Define the headers
set(HEADERS spherical_navier_stokes_elements.h
    refineable_spherical_navier_stokes_elements.h)

# Define the sources
set(SOURCES spherical_navier_stokes_elements.cc
    refineable_spherical_navier_stokes_elements.cc)

# Set the name of the library
set(LIBNAME spherical_navier_stokes)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving spherical_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
