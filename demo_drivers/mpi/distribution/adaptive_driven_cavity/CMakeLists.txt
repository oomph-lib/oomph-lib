# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(adaptive_driven_cavity C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

oomph_add_executable(
  NAME adaptive_driven_cavity_load_balance
  SOURCES adaptive_driven_cavity_load_balance.cc
  LIBRARIES oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.adaptive_driven_cavity
  TARGET_DEPENDENCIES adaptive_driven_cavity adaptive_driven_cavity_load_balance
  EXTRA_REQUIRES adaptive_cavity_1_partition.dat adaptive_cavity_2_partition.dat
                 load_balance.mcr prune_and_balance.mcr
  LABELS "mpi;adaptive_driven_cavity;adapt"
  REQUIRES_MPI_RUN_COMMAND)

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
#
# TODO: TEST TO SEE IF I NEED TO USE:
#
# set_tests_properties(mpi.adaptive_driven_cavity PROPERTIES RUN_SERIAL TRUE)
set_tests_properties(mpi.adaptive_driven_cavity
                     PROPERTIES RESOURCE_LOCK REQUIRES_SOLE_USE_OF_CPU_MPI)

# TODO: Patch
if(APPLE)
  set_tests_properties(mpi.adaptive_driven_cavity PROPERTIES DISABLED YES)
endif()

message(VERBOSE "Leaving adaptive_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
