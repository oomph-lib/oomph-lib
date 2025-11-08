# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered polar_navier_stokes subdirectory")

# Define the headers
set(HEADERS polar_navier_stokes_elements.h
    refineable_polar_navier_stokes_elements.h polar_fluid_traction_elements.h
    polar_stress_integral_elements.h)

# Define the sources
set(SOURCES polar_navier_stokes_elements.cc
    refineable_polar_navier_stokes_elements.cc)

# Set the name of the library
set(LIBNAME polar_navier_stokes)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving polar_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
