# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered poroelasticity subdirectory")

# Define the headers
set(HEADERS elasticity_tensor.h poroelasticity_elements.h
    Tporoelasticity_elements.h poroelasticity_face_elements.h)

# Define the sources
set(SOURCES elasticity_tensor.cc poroelasticity_elements.cc
    Tporoelasticity_elements.cc)

# Set the name of the library
set(LIBNAME poroelasticity)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving poroelasticity subdirectory")
# ------------------------------------------------------------------------------
