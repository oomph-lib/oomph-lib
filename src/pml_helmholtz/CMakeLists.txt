# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pml_helmholtz subdirectory")

# Define the headers
set(HEADERS
    pml_helmholtz_elements.h
    pml_helmholtz_flux_elements.h
    Tpml_helmholtz_elements.h
    refineable_pml_helmholtz_elements.h
    complex_smoother.h
    helmholtz_geometric_multigrid.h)

# Define the sources
set(SOURCES pml_helmholtz_elements.cc refineable_pml_helmholtz_elements.cc)

# Set the name of the library
set(LIBNAME pml_helmholtz)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving pml_helmholtz subdirectory")
# ------------------------------------------------------------------------------
