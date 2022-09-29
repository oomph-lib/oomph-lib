# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered demo_drivers subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(demo_drivers C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

# Add code coverage compiler flags, if requested
if(OOMPH_ENABLE_CODE_COVERAGE AND NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  oomph_enable_code_coverage()
endif()

set(SUBDIRS
    FAQ
    advection_diffusion
    all_foeppl_von_karman
    axisym_advection_diffusion
    axisym_displ_based_foeppl_von_karman
    axisym_foeppl_von_karman
    axisym_linear_elasticity
    axisym_navier_stokes
    axisym_poroelasticity
    beam
    bifurcation_tracking
    biharmonic
    darcy
    eigenproblems
    extrude_quad_mesh_to_cube_mesh
    flux_transport
    foeppl_von_karman
    fourier_decomposed_helmholtz
    generalised_newtonian_axisym_navier_stokes
    generalised_newtonian_navier_stokes
    heat_transfer_and_melting
    helmholtz
    homogenisation
    interaction
    linear_elasticity
    linear_solvers
    linear_wave
    linearised_axisym_navier_stokes
    linearised_free_surface_axisym_navier_stokes
    meshing
    multi_physics
    multigrid
    navier_stokes
    ode
    optimisation
    pml_fourier_decomposed_helmholtz
    pml_helmholtz
    pml_time_harmonic_linear_elasticity
    poisson
    polar_navier_stokes
    reaction_diffusion
    shell
    solid
    spherical_advection_diffusion
    spherical_navier_stokes
    steady_axisym_advection_diffusion
    time_harmonic_fourier_decomposed_linear_elasticity
    time_harmonic_linear_elasticity
    unsteady_heat
    womersley
    young_laplace)

if(OOMPH_HAS_MPI)
  list(APPEND SUBDIRS mpi)
endif()

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving demo_drivers subdirectory")
# ------------------------------------------------------------------------------
