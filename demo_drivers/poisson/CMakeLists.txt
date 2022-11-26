# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered poisson subdirectory")

set(SUBDIRS
    poisson_with_singularity
    one_d_poisson
    one_d_poisson_generic_only
    one_d_poisson_adapt
    one_d_poisson_hp_adapt
    two_d_poisson
    two_d_poisson_flux_bc
    two_d_poisson_flux_bc2
    two_d_poisson_flux_bc_adapt
    fish_poisson
    two_d_poisson_adapt
    two_d_poisson_hp_adapt
    fish_poisson2
    elastic_poisson
    eighth_sphere_poisson
    eighth_sphere_poisson_hp_adapt
    spectral)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving poisson subdirectory")
# ------------------------------------------------------------------------------
