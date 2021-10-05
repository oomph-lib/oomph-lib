# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered hp_adaptive_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(hp_adaptive_driven_cavity C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME hp_adaptive_driven_cavity
  SOURCES hp_adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes MPI::MPI_CXX
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.hp_adaptive_driven_cavity
  TARGET_DEPENDENCIES hp_adaptive_driven_cavity
  EXTRA_REQUIRES hp_adaptive_cavity_partition.dat
  LABELS "mpi;hp_adaptive_driven_cavity;adapt"
  REQUIRES_MPI_RUN_COMMAND)

# Indicate that this test requires the resource REQUIRES_SOLE_USE_OF_CPU_MPI.
# We'll pass this to all other MPI-enabled self-tests to ensure that parallel
# self-test do NOT run at the same time!
set_tests_properties(mpi.hp_adaptive_driven_cavity
                     PROPERTIES RESOURCE_LOCK REQUIRES_SOLE_USE_OF_CPU_MPI)

# TODO: Patch
if(APPLE)
  set_tests_properties(mpi.hp_adaptive_driven_cavity PROPERTIES DISABLED YES)
endif()

message(VERBOSE "Leaving hp_adaptive_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
