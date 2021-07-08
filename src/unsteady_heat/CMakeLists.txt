# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_heat subdirectory")

# Define the headers
set(HEADERS unsteady_heat_elements.h unsteady_heat_flux_elements.h
    refineable_unsteady_heat_elements.h Tunsteady_heat_elements.h)

# Define the sources
set(SOURCES unsteady_heat_elements.cc refineable_unsteady_heat_elements.cc
    Tunsteady_heat_elements.cc)

# Set the name of the library
set(LIBNAME unsteady_heat)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
