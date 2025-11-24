# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pml_time_harmonic_linear_elasticity subdirectory")

# Define the headers
set(HEADERS
    pml_time_harmonic_elasticity_tensor.h
    pml_time_harmonic_linear_elasticity_elements.h
    pml_time_harmonic_linear_elasticity_traction_elements.h
    Tpml_time_harmonic_linear_elasticity_elements.h)

# Define the sources
set(SOURCES pml_time_harmonic_elasticity_tensor.cc
    pml_time_harmonic_linear_elasticity_elements.cc
    Tpml_time_harmonic_linear_elasticity_elements.cc)

# Set the name of the library
set(LIBNAME pml_time_harmonic_linear_elasticity)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving pml_time_harmonic_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
