# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered helmholtz subdirectory")

# Define the headers
set(HEADERS helmholtz_elements.h helmholtz_flux_elements.h
    refineable_helmholtz_elements.h Thelmholtz_elements.h
    helmholtz_bc_elements.h)

# Define the sources
set(SOURCES helmholtz_elements.cc refineable_helmholtz_elements.cc
    Thelmholtz_elements.cc)

# Set the name of the library
set(LIBNAME helmholtz)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving helmholtz subdirectory")
# ------------------------------------------------------------------------------
