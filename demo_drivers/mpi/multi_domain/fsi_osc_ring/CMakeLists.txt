# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_osc_ring subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fsi_osc_ring C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_osc_ring
  SOURCES fsi_osc_ring.cc
  LIBRARIES oomph::beam oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.fsi_osc_ring
  TARGET_DEPENDENCIES fsi_osc_ring
  EXTRA_REQUIRES fsi_osc_ring_partition.dat
  LABELS mpi fsi_osc_ring
  REQUIRES_MPI_RUN_COMMAND)

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
set_tests_properties(mpi.fsi_osc_ring PROPERTIES RESOURCE_LOCK
                                                 REQUIRES_SOLE_USE_OF_CPU_MPI)

message(VERBOSE "Leaving fsi_osc_ring subdirectory")
# ------------------------------------------------------------------------------
