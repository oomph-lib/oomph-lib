# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered airy_cantilever subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(airy_cantilever C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME airy_cantilever2_adapt
  SOURCES airy_cantilever2.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_DEFINITIONS REFINE)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.airy_cantilever
  TARGET_DEPENDENCIES airy_cantilever2_adapt
  EXTRA_REQUIRES
    airy_cantilever_0_partition.dat
    airy_cantilever_1_partition.dat
    airy_cantilever_2_partition.dat
    airy_cantilever_3_partition.dat
    airy_cantilever_4_partition.dat
    airy_cantilever_5_partition.dat
    airy_cantilever_6_partition.dat
    airy_cantilever_7_partition.dat
    airy_cantilever_8_partition.dat
    airy_cantilever_9_partition.dat
  LABELS mpi airy_cantilever
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.airy_cantilever PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving airy_cantilever subdirectory")
# ------------------------------------------------------------------------------
