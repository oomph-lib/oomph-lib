# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_poroelasticity subdirectory")

# Define the headers
set(HEADERS axisym_poroelasticity_elements.h Taxisym_poroelasticity_elements.h
    axisym_poroelasticity_face_elements.h
    axisym_poroelastic_fsi_face_elements.h)

# Define the sources
set(SOURCES axisym_poroelasticity_elements.cc
    Taxisym_poroelasticity_elements.cc)

# Set the name of the library
set(LIBNAME axisym_poroelasticity)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving axisym_poroelasticity subdirectory")
# ------------------------------------------------------------------------------
