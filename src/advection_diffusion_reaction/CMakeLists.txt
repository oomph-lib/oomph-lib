# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered advection_diffusion_reaction subdirectory")

# Define the headers
set(HEADERS advection_diffusion_reaction_elements.h
    refineable_advection_diffusion_reaction_elements.h
    Tadvection_diffusion_reaction_elements.h)
# advection_diffusion_reaction_flux_elements.h
# supg_advection_diffusion_reaction_elements.h

# Define the sources
set(SOURCES advection_diffusion_reaction_elements.cc
    refineable_advection_diffusion_reaction_elements.cc
    Tadvection_diffusion_reaction_elements.cc)

# Set the name of the library
set(LIBNAME advection_diffusion_reaction)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving advection_diffusion_reaction subdirectory")
# ------------------------------------------------------------------------------
