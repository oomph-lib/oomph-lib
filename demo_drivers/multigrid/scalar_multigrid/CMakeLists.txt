# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered scalar_multigrid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(scalar_multigrid C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_tanh_flux_bc
  SOURCES two_d_poisson_tanh_flux_bc.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_poisson_tanh_flux_bc_validate
  SOURCES two_d_poisson_tanh_flux_bc_validate.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unit_cube_poisson
  SOURCES unit_cube_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unit_cube_poisson_validate
  SOURCES unit_cube_poisson_validate.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multigrid.scalar_multigrid
  TARGET_DEPENDENCIES
    two_d_poisson_tanh_flux_bc two_d_poisson_tanh_flux_bc_validate
    unit_cube_poisson unit_cube_poisson_validate
  EXTRA_REQUIRES
  LABELS multigrid scalar_multigrid)

message(VERBOSE "Leaving scalar_multigrid subdirectory")
# ------------------------------------------------------------------------------
