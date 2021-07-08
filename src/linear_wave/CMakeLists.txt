# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linear_wave subdirectory")

# Define the headers
set(HEADERS linear_wave_elements.h linear_wave_flux_elements.h
    refineable_linear_wave_elements.h)

# Define the sources
set(SOURCES linear_wave_elements.cc refineable_linear_wave_elements.cc)

# Set the name of the library
set(LIBNAME linear_wave)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving linear_wave subdirectory")
# ------------------------------------------------------------------------------
