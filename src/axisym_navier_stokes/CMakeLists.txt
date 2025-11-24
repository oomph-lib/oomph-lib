# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_navier_stokes subdirectory")

# Define the headers
set(HEADERS axisym_navier_stokes_elements.h
    refineable_axisym_navier_stokes_elements.h Taxisym_navier_stokes_elements.h
    axisym_fluid_traction_elements.h)

# Define the sources
set(SOURCES axisym_navier_stokes_elements.cc
    refineable_axisym_navier_stokes_elements.cc
    Taxisym_navier_stokes_elements.cc axisym_fluid_traction_elements.cc)

# Set the name of the library
set(LIBNAME axisym_navier_stokes)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
