# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_cantilever subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(three_d_cantilever C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_cantilever_adapt
  SOURCES three_d_cantilever.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_DEFINITIONS REFINE
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.three_d_cantilever
  TARGET_DEPENDENCIES three_d_cantilever_adapt
  EXTRA_REQUIRES
    three_d_cantilever_0_partition.dat
    three_d_cantilever_1_partition.dat
    three_d_cantilever_3_partition.dat
    three_d_cantilever_4_partition.dat
    three_d_cantilever_5_partition.dat
    three_d_cantilever_6_partition.dat
    three_d_cantilever_8_partition.dat
    three_d_cantilever_9_partition.dat
    comp.mcr
  LABELS mpi three_d_cantilever
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.three_d_cantilever PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving three_d_cantilever subdirectory")
# ------------------------------------------------------------------------------
