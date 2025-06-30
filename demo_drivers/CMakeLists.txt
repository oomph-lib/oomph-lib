# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered demo_drivers subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(demo_drivers C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../install")
endif()

include(CTest)

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
    gzip
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
    self_test
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

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving demo_drivers subdirectory")
# ------------------------------------------------------------------------------
