# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered circular_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(circular_driven_cavity C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME circular_driven_cavity
  SOURCES circular_driven_cavity.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.circular_driven_cavity
  TARGET_DEPENDENCIES circular_driven_cavity
  EXTRA_REQUIRES circular_cavity_1_partition.dat circular_cavity_2_partition.dat
  LABELS mpi circular_driven_cavity
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.circular_driven_cavity PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving circular_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
