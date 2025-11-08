# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered multi_physics subdirectory")

# Define the headers
set(HEADERS
    fsi_preconditioners.h
    segregated_fsi_solver.h
    pseudo_elastic_preconditioner.h
    pseudo_elastic_fsi_preconditioner.h
    multi_domain_boussinesq_elements.h
    boussinesq_elements.h
    helmholtz_time_harmonic_linear_elasticity_interaction.h
    fourier_decomposed_helmholtz_time_harmonic_linear_elasticity_interaction.h
    pml_helmholtz_time_harmonic_linear_elasticity_interaction.h)
# old_multi_domain_boussinesq_elements.h

# Define the sources
set(SOURCES segregated_fsi_solver.cc pseudo_elastic_preconditioner.cc
    pseudo_elastic_fsi_preconditioner.cc)

# Set the name of the library
set(LIBNAME multi_physics)

# We depend on the generic library
set(LINKLIBS generic)

# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  LINKLIBS ${LINKLIBS})

message(VERBOSE "Leaving multi_physics subdirectory")
# ------------------------------------------------------------------------------
