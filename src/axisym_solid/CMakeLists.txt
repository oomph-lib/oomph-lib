# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_solid subdirectory")

# Define the headers
set(HEADERS
    axisym_cylindrical_solid_elements.h
    refineable_axisym_cylindrical_solid_elements.h
    axisym_cylindrical_solid_with_pressure_elements.h
    refineable_axisym_cylindrical_solid_with_pressure_elements.h
    axisym_cylindrical_solid_traction_elements.h)

# Define the sources
set(SOURCES axisym_cylindrical_solid_elements.cc
    refineable_axisym_cylindrical_solid_elements.cc)

# Set the name of the library
set(LIBNAME axisym_solid)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving axisym_solid subdirectory")
# ------------------------------------------------------------------------------
