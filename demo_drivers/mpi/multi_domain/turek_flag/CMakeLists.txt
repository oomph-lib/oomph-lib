# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered turek_flag subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(turek_flag C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME turek_flag
  SOURCES turek_flag.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME turek_flag_load_balance
  SOURCES turek_flag_load_balance.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.turek_flag
  TARGET_DEPENDENCIES turek_flag turek_flag_load_balance
  EXTRA_REQUIRES turek_flag_partition.dat load_balance.mcr
  LABELS mpi turek_flag
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.turek_flag PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving turek_flag subdirectory")
# ------------------------------------------------------------------------------
