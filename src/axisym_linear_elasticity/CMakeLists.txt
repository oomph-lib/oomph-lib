# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_linear_elasticity subdirectory")

# Define the headers
set(HEADERS axisym_linear_elasticity_elements.h
    axisym_linear_elasticity_traction_elements.h
    Taxisym_linear_elasticity_elements.h)

# Define the sources
set(SOURCES axisym_linear_elasticity_elements.cc
    Taxisym_linear_elasticity_elements.cc)

# Set the name of the library
set(LIBNAME axisym_linear_elasticity)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving axisym_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
