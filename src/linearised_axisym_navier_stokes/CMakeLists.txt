# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linearised_axisym_navier_stokes subdirectory")

# Define the headers
set(HEADERS linearised_axisym_navier_stokes_elements.h
    refineable_linearised_axisym_navier_stokes_elements.h
    multi_domain_linearised_axisym_navier_stokes_elements.h)

# Define the sources
set(SOURCES linearised_axisym_navier_stokes_elements.cc
    refineable_linearised_axisym_navier_stokes_elements.cc)

# Set the name of the library
set(LIBNAME linearised_axisym_navier_stokes)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving linearised_axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
