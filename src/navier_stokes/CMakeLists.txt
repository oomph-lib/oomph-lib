# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered navier_stokes subdirectory")

# Define the headers
set(HEADERS
    navier_stokes_elements.h
    fluid_traction_elements.h
    refineable_navier_stokes_elements.h
    Tnavier_stokes_elements.h
    navier_stokes_preconditioners.h
    navier_stokes_surface_power_elements.h
    navier_stokes_surface_drag_torque_elements.h
    navier_stokes_flux_control_elements.h
    impose_parallel_outflow_element.h
    impose_impenetrability_element.h
    lagrange_enforced_flow_preconditioner.h
    vorticity_smoother.h)

# Define the sources
set(SOURCES navier_stokes_elements.cc refineable_navier_stokes_elements.cc
    Tnavier_stokes_elements.cc navier_stokes_preconditioners.cc
    lagrange_enforced_flow_preconditioner.cc)

# Set the name of the library
set(LIBNAME navier_stokes)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving navier_stokes subdirectory")
# ------------------------------------------------------------------------------
