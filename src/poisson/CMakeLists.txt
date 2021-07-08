# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered poisson subdirectory")

# Define the headers
set(HEADERS
    poisson_elements.h
    poisson_flux_elements.h
    spectral_poisson_elements.h
    refineable_poisson_elements.h
    refineable_spectral_poisson_elements.h
    Tpoisson_elements.h)

# Define the sources
set(SOURCES poisson_elements.cc spectral_poisson_elements.cc
    refineable_poisson_elements.cc Tpoisson_elements.cc)

# Set the name of the library
set(LIBNAME poisson)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving poisson subdirectory")
# ------------------------------------------------------------------------------
