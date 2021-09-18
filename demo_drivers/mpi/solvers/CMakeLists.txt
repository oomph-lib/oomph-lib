# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered solvers subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(solvers C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME driven_cavity
  SOURCES driven_cavity.cc
  LIBRARIES oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

oomph_add_executable(
  NAME direct_solver_test
  SOURCES direct_solver_test.cc
  LIBRARIES oomph::poisson MPI::MPI_CXX
  CXX_STANDARD 14)

oomph_add_executable(
  NAME airy_cantilever
  SOURCES airy_cantilever.cc
  LIBRARIES oomph::solid oomph::constitutive MPI::MPI_CXX
  CXX_STANDARD 14)

oomph_add_executable(
  NAME fsi_channel_with_leaflet
  SOURCES fsi_channel_with_leaflet.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::generic MPI::MPI_CXX
  CXX_STANDARD 14)

if(OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME TrilinosSolver_test
    SOURCES TrilinosSolver_test.cc
    LIBRARIES oomph::poisson MPI::MPI_CXX
    CXX_STANDARD 14)

  oomph_add_executable(
    NAME two_d_multi_poisson
    SOURCES two_d_multi_poisson.cc
    LIBRARIES oomph::generic MPI::MPI_CXX
    CXX_STANDARD 14)
endif()

# ---------------------------------- TESTING ----------------------------------

set(TEST_TARGET_DEPS driven_cavity direct_solver_test airy_cantilever
    fsi_channel_with_leaflet)

if(OOMPH_HAS_TRILINOS)
  list(APPEND TEST_TARGET_DEPS TrilinosSolver_test two_d_multi_poisson)
endif()

oomph_add_test(
  TEST_NAME mpi.solvers
  TARGET_DEPENDENCIES ${TEST_TARGET_DEPS}
  EXTRA_REQUIRES multi_poisson_partition.dat
  LABELS "mpi;solvers"
  REQUIRES_MPI_RUN_COMMAND)

message(VERBOSE "Leaving solvers subdirectory")
# ------------------------------------------------------------------------------
