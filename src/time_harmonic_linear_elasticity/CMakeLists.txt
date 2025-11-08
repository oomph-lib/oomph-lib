# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered time_harmonic_linear_elasticity subdirectory")

# Define the headers
set(HEADERS
    time_harmonic_elasticity_tensor.h
    time_harmonic_linear_elasticity_elements.h
    time_harmonic_linear_elasticity_traction_elements.h
    Ttime_harmonic_linear_elasticity_elements.h
    refineable_time_harmonic_linear_elasticity_elements.h)

# Define the sources
set(SOURCES
    time_harmonic_elasticity_tensor.cc
    time_harmonic_linear_elasticity_elements.cc
    refineable_time_harmonic_linear_elasticity_elements.cc
    Ttime_harmonic_linear_elasticity_elements.cc)

# Set the name of the library
set(LIBNAME time_harmonic_linear_elasticity)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving time_harmonic_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
