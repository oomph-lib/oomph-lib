# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linear_elasticity subdirectory")

# Define the headers
set(HEADERS elasticity_tensor.h linear_elasticity_elements.h
    linear_elasticity_traction_elements.h Tlinear_elasticity_elements.h
    refineable_linear_elasticity_elements.h)

# Define the sources
set(SOURCES elasticity_tensor.cc linear_elasticity_elements.cc
    refineable_linear_elasticity_elements.cc Tlinear_elasticity_elements.cc
    linear_elasticity_traction_elements.cc)

# Set the name of the library
set(LIBNAME linear_elasticity)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
