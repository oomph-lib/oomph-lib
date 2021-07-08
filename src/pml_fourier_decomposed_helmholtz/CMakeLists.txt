# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pml_fourier_decomposed_helmholtz subdirectory")

# Define the headers
set(HEADERS pml_fourier_decomposed_helmholtz_elements.h
    pml_fourier_decomposed_helmholtz_flux_elements.h
    Tpml_fourier_decomposed_helmholtz_elements.h)

# Define the sources
set(SOURCES pml_fourier_decomposed_helmholtz_elements.cc
    Tpml_fourier_decomposed_helmholtz_elements.cc)

# Set the name of the library
set(LIBNAME pml_fourier_decomposed_helmholtz)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving pml_fourier_decomposed_helmholtz subdirectory")
# ------------------------------------------------------------------------------
