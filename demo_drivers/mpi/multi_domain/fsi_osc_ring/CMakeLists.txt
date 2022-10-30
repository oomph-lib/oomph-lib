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
  LIBRARIES oomph::beam oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.fsi_osc_ring
  TARGET_DEPENDENCIES fsi_osc_ring
  EXTRA_REQUIRES fsi_osc_ring_partition.dat
  LABELS mpi fsi_osc_ring
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.fsi_osc_ring PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving fsi_osc_ring subdirectory")
# ------------------------------------------------------------------------------
