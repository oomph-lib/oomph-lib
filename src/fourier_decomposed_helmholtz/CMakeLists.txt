# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fourier_decomposed_helmholtz subdirectory")

# Define the headers
set(HEADERS
    fourier_decomposed_helmholtz_elements.h
    fourier_decomposed_helmholtz_flux_elements.h
    fourier_decomposed_helmholtz_bc_elements.h
    Tfourier_decomposed_helmholtz_elements.h)

# Define the sources
set(SOURCES fourier_decomposed_helmholtz_elements.cc
    Tfourier_decomposed_helmholtz_elements.cc)

# Set the name of the library
set(LIBNAME fourier_decomposed_helmholtz)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving fourier_decomposed_helmholtz subdirectory")
# ------------------------------------------------------------------------------
