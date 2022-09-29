# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_poisson_flux_bc_adapt subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_poisson_flux_bc_adapt C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_flux_bc_adapt
  SOURCES two_d_poisson_flux_bc_adapt.cc
  LIBRARIES oomph::poisson MPI::MPI_CXX
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.two_d_poisson_flux_bc_adapt
  TARGET_DEPENDENCIES two_d_poisson_flux_bc_adapt
  EXTRA_REQUIRES two_d_poisson_flux_partition.dat
  LABELS "mpi;two_d_poisson_flux_bc_adapt;adapt"
  REQUIRES_MPI_RUN_COMMAND)

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
set_tests_properties(mpi.two_d_poisson_flux_bc_adapt
                     PROPERTIES RESOURCE_LOCK REQUIRES_SOLE_USE_OF_CPU_MPI)

message(VERBOSE "Leaving two_d_poisson_flux_bc_adapt subdirectory")
# ------------------------------------------------------------------------------
