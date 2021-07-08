# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered solid subdirectory")

# Define the headers
set(HEADERS solid_elements.h refineable_solid_elements.h
    solid_traction_elements.h solid_preconditioners.h)

# Define the sources
set(SOURCES solid_elements.cc refineable_solid_elements.cc
    solid_preconditioners.cc)

# Set the name of the library
set(LIBNAME solid)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving solid subdirectory")
# ------------------------------------------------------------------------------
