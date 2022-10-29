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
  LIBRARIES oomph::solid oomph::constitutive MPI::MPI_CXX
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

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
set_tests_properties(mpi.three_d_cantilever
                     PROPERTIES RESOURCE_LOCK REQUIRES_SOLE_USE_OF_CPU_MPI)

# TODO: Patch
if(APPLE)
  set_tests_properties(mpi.three_d_cantilever PROPERTIES DISABLED YES)
endif()

message(VERBOSE "Leaving three_d_cantilever subdirectory")
# ------------------------------------------------------------------------------
