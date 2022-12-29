# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mpi subdirectory")

set(SUBDIRS
    general_mpi
    adaptive_driven_cavity
    adaptive_driven_cavity_load_balance
    boussinesq_convection
    fsi_channel_with_leaflet
    turek_flag
    two_d_poisson_flux_bc_adapt
    block_preconditioners
    distributed_general_purpose_block_preconditioners
    distributed_linear_algebra_infrastructure)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving mpi subdirectory")
# ------------------------------------------------------------------------------
