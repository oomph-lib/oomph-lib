# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered steady_axisym_advection_diffusion subdirectory")

# Define the headers
set(HEADERS steady_axisym_advection_diffusion_elements.h)

# Define the sources
set(SOURCES steady_axisym_advection_diffusion_elements.cc)

# Set the name of the library
set(LIBNAME steady_axisym_advection_diffusion)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving steady_axisym_advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
