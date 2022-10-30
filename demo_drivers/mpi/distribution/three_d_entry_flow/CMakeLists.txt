# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_entry_flow subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(three_d_entry_flow C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_entry_flow
  SOURCES three_d_entry_flow.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.three_d_entry_flow
  TARGET_DEPENDENCIES three_d_entry_flow
  EXTRA_REQUIRES three_d_entry_flow_1_partition.dat
                 three_d_entry_flow_2_partition.dat
  LABELS mpi three_d_entry_flow
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.three_d_entry_flow PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving three_d_entry_flow subdirectory")
# ------------------------------------------------------------------------------
