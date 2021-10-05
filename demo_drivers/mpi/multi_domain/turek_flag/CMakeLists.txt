# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered turek_flag subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(turek_flag C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME turek_flag
  SOURCES turek_flag.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

oomph_add_executable(
  NAME turek_flag_load_balance
  SOURCES turek_flag_load_balance.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.turek_flag
  TARGET_DEPENDENCIES turek_flag turek_flag_load_balance
  EXTRA_REQUIRES turek_flag_partition.dat load_balance.mcr
  LABELS "mpi;turek_flag"
  REQUIRES_MPI_RUN_COMMAND)

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
set_tests_properties(mpi.turek_flag PROPERTIES RESOURCE_LOCK
                                               REQUIRES_SOLE_USE_OF_CPU_MPI)

# TODO: Patch
if(APPLE)
  set_tests_properties(mpi.turek_flag PROPERTIES DISABLED YES)
endif()

message(VERBOSE "Leaving turek_flag subdirectory")
# ------------------------------------------------------------------------------
