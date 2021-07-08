# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered generalised_newtonian_navier_stokes subdirectory")

# Define the headers
set(HEADERS generalised_newtonian_navier_stokes_elements.h
    generalised_newtonian_Tnavier_stokes_elements.h
    generalised_newtonian_refineable_navier_stokes_elements.h)

# Define the sources
set(SOURCES generalised_newtonian_navier_stokes_elements.cc
    generalised_newtonian_Tnavier_stokes_elements.cc
    generalised_newtonian_refineable_navier_stokes_elements.cc)

# Set the name of the library
set(LIBNAME generalised_newtonian_navier_stokes)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving generalised_newtonian_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
