# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered hp_adaptive_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(hp_adaptive_driven_cavity C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME hp_adaptive_driven_cavity
  SOURCES hp_adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.hp_adaptive_driven_cavity
  TARGET_DEPENDENCIES hp_adaptive_driven_cavity
  EXTRA_REQUIRES hp_adaptive_cavity_partition.dat
  LABELS mpi hp_adaptive_driven_cavity adapt
  REQUIRES_MPI_RUN_COMMAND)

# Don't run this test alongside any other tests
set_tests_properties(mpi.hp_adaptive_driven_cavity PROPERTIES RUN_SERIAL TRUE)

message(VERBOSE "Leaving hp_adaptive_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
