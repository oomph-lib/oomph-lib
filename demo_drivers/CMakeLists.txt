# cmake-format: off
# ==============================================================================
# Example driver codes for testing/illustrating oomph-lib functionality. Go wild.
# ==============================================================================
# cmake-format: on
cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)

# # TODO: Find out how else to fix the issue with deployment targets causing a
# lot # of warnings on OSX. if(APPLE) set(CMAKE_OSX_DEPLOYMENT_TARGET "11.5"
# CACHE INTERNAL "") endif()

project(demo_drivers C CXX Fortran)
find_package(oomphlib REQUIRED)
enable_testing()
include(CTest)

# Helpful info. for the user
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered demo_drivers subdirectory")

set(SUBDIRS
    poisson
    unsteady_heat
    ode
    heat_transfer_and_melting
    advection_diffusion
    helmholtz
    time_harmonic_linear_elasticity
    axisym_advection_diffusion
    spherical_advection_diffusion
    steady_axisym_advection_diffusion
    navier_stokes
    axisym_navier_stokes
    polar_navier_stokes
    spherical_navier_stokes
    solid
    beam
    shell
    linear_wave
    eigenproblems
    interaction
    meshing
    multi_physics
    optimisation
    bifurcation_tracking
    FAQ
    linear_solvers
    linear_elasticity
    womersley
    biharmonic
    reaction_diffusion
    flux_transport
    young_laplace
    time_harmonic_fourier_decomposed_linear_elasticity
    fourier_decomposed_helmholtz
    foeppl_von_karman
    axisym_foeppl_von_karman
    axisym_displ_based_foeppl_von_karman
    pml_helmholtz
    axisym_linear_elasticity
    pml_fourier_decomposed_helmholtz
    homogenisation
    darcy
    axisym_poroelasticity
    pml_time_harmonic_linear_elasticity
    linearised_axisym_navier_stokes
    linearised_free_surface_axisym_navier_stokes
    generalised_newtonian_navier_stokes
    generalised_newtonian_axisym_navier_stokes
    multigrid
    all_foeppl_von_karman
    extrude_quad_mesh_to_cube_mesh)

if(OOMPH_HAS_MPI)
  list(APPEND SUBDIRS mpi)
endif()

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()
# add_subdirectory(poisson)

message(VERBOSE "Leaving demo_drivers subdirectory")
