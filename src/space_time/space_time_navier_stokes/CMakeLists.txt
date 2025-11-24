# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered space_time_navier_stokes subdirectory")

# Define the headers
set(HEADERS
    mixed_order_petrov_galerkin_space_time_fluid_traction_elements.h
    mixed_order_petrov_galerkin_space_time_navier_stokes_elements.h
    refineable_mixed_order_petrov_galerkin_space_time_navier_stokes_elements.h
    mixed_order_petrov_galerkin_space_time_block_preconditionable_elements.h)

# Define the sources
set(SOURCES mixed_order_petrov_galerkin_space_time_navier_stokes_elements.cc
    refineable_mixed_order_petrov_galerkin_space_time_navier_stokes_elements.cc
    mixed_order_petrov_galerkin_space_time_block_preconditionable_elements.cc)

# Set the name of the library
set(LIBNAME space_time_navier_stokes)

# We depend on the generic library
set(LINKLIBS generic space_time_block_preconditioner)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving space_time_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
